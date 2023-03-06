# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 09/07/2022 22:58
@Author: yao
"""

import os
import argparse

import re

from collections import defaultdict


def main(bioc_file: str, save_file: str):

    process_count = 0
    null = 'null'
    with open(bioc_file) as f, open(save_file, 'w') as wf:
        wf.write('PMID\tPMCID\t'
                 'Title\tJournal\t'
                 'Year\n')
        for line in f:

            process_count += 1
            if process_count % 1000 == 0:
                print(f'Processed: {process_count:,}')

            lit_dic = eval(line.strip())

            # for _key in lit_dic.keys():
            #     print(_key)
            #     print(lit_dic[_key])
            #     input()

            pmid = lit_dic['pmid']
            pmcid = lit_dic['pmcid']

            title = lit_dic['passages'][0]['text']

            # print(lit_dic['passages'])

            journal = lit_dic['journal']
            year = lit_dic['year']

            if re.findall(r'\d+', journal):
                # print(journal)
                num = re.findall(r'\d+', journal)[0]
                journal = journal[:journal.find(num)].strip().strip('.')

            wf.write(f'{pmid}\t{pmcid}\t'
                     f'{title}\t{journal}\t{year}\n')
    print(f'{save_file} save done.')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Get_LitInfo_from_Bioc')
    parser.add_argument('-if', dest='bioc_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    main(args.bioc_file, args.save_file)



