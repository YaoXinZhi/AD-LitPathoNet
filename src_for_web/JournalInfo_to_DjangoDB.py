# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 04/08/2022 19:43
@Author: yao
"""

import argparse

from collections import defaultdict

class JournalInFo:
    def __init__(self, ISSN: str, JournalName: str, JournalAbbr: str, ):
        self.ISSN = ISSN
        self.JournalName = JournalName
        self.JournalAbbr = JournalAbbr

        self.ImpactFactor = ''
        self.HIndex = ''
        self.Altmetric = ''


def read_journal_info_file(journal_file: str):
    pmid_to_pub_year = {}
    pmid_to_issn = {}
    issn_to_journal = defaultdict(JournalInFo)

    with open(journal_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')

            if len(l) < 8:
                continue

            pmid = l[ 0 ]
            journal_name = l[ 4 ]
            journal_abbr = l[ 5 ]
            issn = l[ 6 ]

            pub_data = l[7]

            pub_year = pub_data.split('-')[ 0 ]
            pmid_to_pub_year[ pmid ] = pub_year

            pmid_to_issn[ pmid ] = issn
            issn_to_journal[ issn ] = JournalInFo(issn, journal_name, journal_abbr)

    print(f'Total PMID: {len(pmid_to_issn):,}.')
    print(f'Total ISSN: {len(issn_to_journal):,}.')
    return pmid_to_pub_year, pmid_to_issn, issn_to_journal


def read_impact_factor_file(impact_factor_file: str, issn_to_journal: dict):
    updated_issn_set = set()
    with open(impact_factor_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            if len(l)< 10:
                continue

            issn = l[ 6 ]
            impact_factor = l[ 8 ]
            h_index = l[ 9 ]

            issn_to_journal[ issn ].ImpactFactor = impact_factor
            issn_to_journal[ issn ].HIndex = h_index

            updated_issn_set.add(issn)

    print(f'{len(updated_issn_set):,} journals have ImpactFactor/H-index.')
    return issn_to_journal


def read_alt_metric_file(alt_metric_file: str, issn_to_journal: dict):
    updated_issn_set = set()
    with open(alt_metric_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')

            if len(l) < 9:
                continue

            issn = l[ 6 ]
            alt_metric = l[8]

            issn_to_journal[ issn ].Altmetric = alt_metric
            updated_issn_set.add(issn)

    print(f'{len(updated_issn_set):,} journals have Altmetric.')
    return issn_to_journal


def main(journal_info_file: str,
         impact_factor_file: str, alt_metric_file: str,
         save_path: str, prefix: str):
    """
    生成两个文件
    1. 每条证据的发表年份, ISSN （pmid to PublicationYear/ISSN）
    2. 每个杂志的信息 （去重）
    """

    literature_save_file = f'{save_path}/{prefix}.literature.tsv'
    journal_save_file = f'{save_path}/{prefix}.Journal.tsv'

    pmid_to_pub_year, pmid_to_issn, issn_to_journal = read_journal_info_file(journal_info_file)

    issn_to_journal = read_impact_factor_file(impact_factor_file, issn_to_journal)

    issn_to_journal = read_alt_metric_file(alt_metric_file, issn_to_journal)

    with open(literature_save_file, 'w') as wf:
        wf.write(f'PMID\tISSN\tPublishedYear\n')
        for pmid in pmid_to_pub_year.keys():
            pub_year = pmid_to_pub_year[ pmid ]
            issn = pmid_to_issn[ pmid ]

            wf.write(f'{pmid}\t{issn}\t{pub_year}\n')
    print(f'{literature_save_file} save done.')

    with open(journal_save_file, 'w') as wf:
        wf.write(f'ISSN\tJournalName\tJournalAbbr\t'
                 f'ImpactFactor\tH-index\t'
                 f'Altmetric\n')
        for issn, journal in issn_to_journal.items():
            wf.write(f'{issn}\t{journal.JournalName}\t{journal.JournalAbbr}\t'
                     f'{journal.ImpactFactor}\t{journal.HIndex}\t'
                     f'{journal.Altmetric}\n')
    print(f'{journal_save_file} save done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='JournalInfo_to_DjangoDB.')
    parser.add_argument('-jf', dest='journal_info_file', required=True)
    parser.add_argument('-if', dest='impact_factor_file', required=True)
    parser.add_argument('-af', dest='alt_metric_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    parser.add_argument('-pr', dest='prefix', required=True)
    args = parser.parse_args()

    main(args.journal_info_file, args.impact_factor_file, args.alt_metric_file,
         args.save_path, args.prefix)
