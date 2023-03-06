# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 18/07/2022 12:03
@Author: XINZHI YAO
"""

import re
import time

import argparse
import requests
from tqdm import tqdm

from bs4 import BeautifulSoup

def url_parser(url: str):

    headers = {'User-Agent':'user-agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/80.0.3987.87 Safari/537.36', 'Cookie':'ncbi_sid=47EF8993E4BE3193_0587SID; _ga=GA1.2.780797648.1582031800; pmc.article.report=; books.article.report=; _gid=GA1.2.288894097.1585885708; QSI_HistorySession=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2F~1585922167314; WebEnv=1zmWDqnmAX_v2YHx4rlkrq3Yhg3ZCIGss-UrAFg-SZPikri1ywXU5zjN2MwUgcrtonBWnvUs1EJyzzVzB2wC14pFkl6eh-708Vl%4047EF8993E4BE3193_0587SID; ncbi_pinger=N4IgDgTgpgbg+mAFgSwCYgFwgCIAYBC2ArAMIDMRAHGdmQCwCcAbEbm7mbgwOxEBiDWvzoA6AIwiAtnEogANCACuAOwA2AewCGqZVAAeAF0ygATJhABzAI6KoEAJ7yQZc5tWqndcwGcomiADGiE5E5hCKqlAAvIj2qBCaAD6oAEZRYIopklCoid7qioHRJMgGmgbI6soAylDK+RB5galRAAqZALI5iQb2YNH+FQGRiU5iXli+/kFOJrjmeISkFNS0jCzsHFy8AkJ8ohLSsgomYubWtg4YbqoYU4GIGOGRMXEJyWkZWd35hQHFpXKlRqdQaTQCLXaKS6uV6/Sig2QwygoxOZiwAHcsSJlAEUsgcapJDjkIgRBZ1DBZgxzGIGHMnJxabhKPMFGQzlg6Qz2eiQGJcEQXOyXFgDOEoIyJiBxbZGbIsIyaVgiELBTSFHR5lgmDwGJRQprRSBcCITJQRMKQHRpSoNNpdIZPKEsGyQELaawQtKGFqQkxaUwOCFuLSvAomJz+dxuLIAL5xoA'}

    html = requests.get(url, headers=headers)
    html = html.content.decode('utf-8')
    # soup = BeautifulSoup(html, 'html.parser')
    return html

def read_issn_file(issn_file:str):

    issn_set = set()
    with open(issn_file) as f:
        for line in f:
            l = line.strip().split('\t')
            if len(l)< 3:
                continue
            issn = l[6]
            issn_set.add(issn)
    return issn_set

def get_if(issn_file: str, save_file: str):

    # issn_set = read_issn_file(issn_file)

    err_count = 0
    process_count = 0
    with open(issn_file) as f, open(save_file, 'w') as wf:
        head = f.readline().strip()
        wf.write(f'{head}\tImpactFactor\tH-index\n')
        for line in f:
            l = line.strip().split('\t')

            if len(l) < 3:
                continue
            process_count += 1

            if process_count % 50 == 0:
                time.sleep(10)
                print(f'{process_count:,} ISSN processed.')
            try:
                pmid, pmcid, title, doi, journal, journal_abbr, issn, pub_date = l

                url = f'https://www.medsci.cn/sci/index.do?fullname={issn}&page=1'


                doc = url_parser(url)
            except:
                err_count += 1
                wf.write(f'{line.strip()}\tNone\tNone\n')
                continue

            impact_factor = re.findall(r'<span>影响指数：(\d+\.\d+)</span>', doc)
            if impact_factor:
                impact_factor = impact_factor[0]
            else:
                impact_factor = 'None'

            H_index = re.findall(r'<span>H指数：(\d+)</span>', doc)
            if H_index:
                H_index = H_index[0]
            else:
                H_index = 'None'

            wf.write(f'{line.strip()}\t{impact_factor}\t{H_index}\n')
            # wf.write(f'{issn}\t{impact_factor}\t{H_index}\n')

    print(f'Error: {err_count:,}.')
    print(f'{save_file} save done.')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get ImpactFact From MedSci.')
    parser.add_argument('-jf', dest='JournalInfo_file', required=True)
    parser.add_argument('-sf', dest='IF_save_file', required=True)
    args = parser.parse_args()

    # JournalInfo_file = '../result/ad-total-20220427/LitInfo_dir/ad.JournalInfo.tsv'
    # IF_save_file = '../result/ad-total-20220427/LitInfo_dir/ad.ImpactFactor.MedSci.tsv'

    get_if(args.JournalInfo_file, args.IF_save_file)
