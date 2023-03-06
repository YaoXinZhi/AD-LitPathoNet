# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 22/05/2022 18:55
@Author: XINZHI YAO
"""

import os
import argparse

from tqdm import tqdm

import re
import requests

from easy_entrez import EntrezAPI

from xml.dom import minidom
from xml.etree import ElementTree

def xml_to_string(element):
    return (
        minidom.parseString(ElementTree.tostring(element))
        .toprettyxml(indent=' ' * 4)
    )

def read_rs_file(rs_file: str):

    rs_set = set()
    with open(rs_file) as f:
        for line in f:
            l = line.strip()
            rs_set.add(l)

    return rs_set

def get_entrez(rs_id: str):

    # non-coding mutation: rs60320384
    base_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id={rs_id}&retmode=xml'
    try:
        html = requests.get(base_url)
        doc = html.content
        if type(doc) == bytes:
            doc = bytes.decode(doc)
         
        #print(doc)
        #print(re.findall(r'<SYMBOL>(.*?)</SYMBOL>', doc))
        #input()
        symbol = re.findall(r'<NAME>(.*?)</NAME>', doc)

        tax_id = re.findall(r'<TAX_ID>(.*?)</TAX_ID>', doc)

        entrez = re.findall(r'<GENE_ID>(.*?)</GENE_ID>', doc)

        fxn_class = re.findall(r'<FXN_CLASS>(.*?)</FXN_CLASS>', doc)
        return tax_id, entrez, symbol, fxn_class
    except:
        return '', '', '', ''

def main_for_easy_entrez(rs_file: str, save_file: str):

    rs_set = read_rs_file(rs_file)

    entrez_api = EntrezAPI(
        'your-tool-name',
        'e@mail.com',
        # optional
        return_type='json'
    )
    namespaces = {'ns0': 'https://www.ncbi.nlm.nih.gov/SNP/docsum'}

    # print(xml_to_string(rs6311))
    match_count = 0
    with open(save_file, 'w') as wf:
        wf.write('tax_id\tEntrez\tSymbol\tFXNClass\n')
        for rs_id in tqdm(rs_set):
            try:
                rs_result = entrez_api.fetch([rs_id], max_results=1, database='snp').data[0]

                tax_id = [name.text
                          for name in rs_result.findall('.//ns0:TAX_ID', namespaces)][0]
                entrez = [name.text
                          for name in rs_result.findall('.//ns0:GENE_E/ns0:GENE_ID', namespaces)][0]
                symbol = [name.text
                          for name in rs_result.findall('.//ns0:GENE_E/ns0:NAME', namespaces)][0]

                fxn_class = [name.text
                             for name in rs_result.findall('.//ns0:FXN_CLASS', namespaces)][0]
                #print(f'{rs_id}\t{tax_id}\t{entrez}\t{symbol}\n')
            except:
                continue
            # print(f'{rs_id}\t{entrez}\t{symbol}\t{fxn_class}')
            if tax_id == '9606' and entrez:
                match_count += 1
                wf.write(f'{rs_id}\t{entrez}\t{symbol}\t{fxn_class}\n')
    print(f'{save_file} save done, matched: {match_count:,}/{len(rs_set):,}.')


def main(rs_file: str, save_file: str):

    rs_set = read_rs_file(rs_file)

    match_count = 0
    with open(save_file, 'w') as wf:
        wf.write('rs_id\tEntrez\n')
        for rs_id in tqdm(rs_set):
            tax_id, entrez, symbol, fxn_class = get_entrez(rs_id)

            if tax_id == '9606' and entrez:
                wf.write(f'{rs_id}\t{entrez}\t{symbol}\t{fxn_class}\n')
                match_count += 1

    print(f'{save_file} save done, match rs_id: {match_count}/{len(rs_set)}.')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get entrez from rsID.')

    parser.add_argument('-rf', dest='rs_file',
                        required=True)
    parser.add_argument('-sf', dest='save_file',
                        required=True)
    args = parser.parse_args()

    # main(args.rs_file, args.save_file)
    main_for_easy_entrez(args.rs_file, args.save_file)


