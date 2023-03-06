# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 17/07/2022 16:24
@Author: XINZHI YAO
"""

import os
import argparse

import time
from tqdm import tqdm

from pubmed_mapper import Article


def read_pmid_file(pmid_file: str):
    """
    PMID PMCID
    """
    pmid_set = set()
    with open(pmid_file) as f:
        # f.readline()
        for line in f:
            PMID = line.strip().split('\t')[0]
            pmid_set.add(PMID)

    print(f'PMID: {len(pmid_set):,}')
    return pmid_set

def pmid_to_journal_info(pmid_file:str, save_file: str):

    pmid_set  =read_pmid_file(pmid_file)

    request_count = 0
    error_pmid_set = set()
    with open(save_file, 'w') as wf:
        wf.write('PMID\tPMCID\t'
                 'Title\tdoi\t'
                 'Journal\tJournalAbbr\t'
                 'ISSN\t'
                 'PubDate\n')
        for pmid in tqdm(pmid_set):

            request_count += 1
            if request_count%50 ==0:
                time.sleep(10)

            try:
                article = Article.parse_pmid(pmid)

                doi = 'None'
                pmc_id = 'null'
                for _id in article.ids:
                    if _id.id_type == 'doi':
                        doi = _id.id_value
                    if _id.id_type == 'pmc':
                        pmc_id = _id.id_value

                title = article.title

                journal_name = article.journal.title
                journal_abbr = article.journal.abbr

                issn = article.journal.issn
                pub_date = article.pubdate

                wf.write(f'{pmid}\t{pmc_id}\t'
                         f'{title}\t{doi}\t'
                         f'{journal_name}\t{journal_abbr}\t'
                         f'{issn}\t{pub_date}\n')
            except:
                wf.write(f'{pmid}\tnull\n')
                continue

    print(f'err: {len(error_pmid_set)}')
    print(f'{save_file} save done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get Journal Info from PMID.')
    parser.add_argument('-pf', dest='pmid_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    pmid_to_journal_info(args.pmid_file, args.save_file)

