# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 18/07/2022 16:17
@Author: XINZHI YAO
"""


import os

import argparse

from pyaltmetric import Altmetric


def main(journal_info_file: str, save_file: str):

    alt_metric_parser = Altmetric()

    process_count = 0
    with open(journal_info_file) as f, open(save_file, 'w') as wf:
        head = f.readline().strip()
        wf.write(f'{head}\tAltMetric\n')
        for line in f:
            l = line.strip().split('\t')

            process_count += 1
            if process_count % 100 == 0:
                print(f'{process_count:,} pmid processed.')

            try:
                pmid, pmcid, title, doi, journal, journal_abbr, issn, pub_date = l

                alt_metric_result = alt_metric_parser.pmid(pmid)

                if alt_metric_result:
                    alt_metric_score = alt_metric_result['score']
                else:
                    alt_metric_score = 'None'

                wf.write(f'{line.strip()}\t{alt_metric_score}\n')
            except:
                wf.write(line)
                print(f'Error: {pmid}')
                continue

    print(f'{save_file} save done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get AltMetric form pmid.')
    parser.add_argument('-jf', dest='journal_info_file', required=True)

    parser.add_argument('-sf', dest='save_file', required=True)

    args = parser.parse_args()

    main(args.journal_info_file, args.save_file)



