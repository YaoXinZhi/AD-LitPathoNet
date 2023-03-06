# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 01/11/2022 10:35
@Author: XINZHI YAO
"""
import os.path

"""
该代码用于筛选原始的chain文件
从而用于生成PNRLE的base-DB-file

并且统计没有标准化的基因突变和概念
"""

import argparse
from collections import defaultdict

class chain_selector:

    def __init__(self):
        self.no_norm_gene = defaultdict(int)
        self.no_norm_var = defaultdict(int)
        self.no_norm_concept = defaultdict(int)

    # @staticmethod
    def chain_filter(self, chain: str):

        # 基本元素都要有
        if len(chain) < 7:
            # Gene ThemeOf Var CauseOf Reg ThemeOf BP
            return False

        # 前三个必须是 Gene ThemeOf Var
        if eval(chain[ 0 ])[ 1 ] != 'Gene' \
                or chain[ 1 ] != 'ThemeOf' \
                or eval(chain[ 2 ])[ 1 ] != 'Var':
            return False

        for element in chain:
            if element in {'CauseOf', 'ThemeOf'}:
                continue
            element = eval(element)
            # 非调控元素必须标准化
            if element[ 1 ] not in {'Reg', 'PosReg', 'NegReg'} and element[ -1 ] == '-':
                if element[ 1 ] == 'Gene':
                    self.no_norm_gene[ element[ 0 ] ] += 1
                elif element[ 1 ] == 'Var':
                    self.no_norm_var[ element[ 0 ] ] += 1
                else:
                    self.no_norm_concept[ element[ 0 ] ] += 1
                return False

        return True

    def main(self, chain_file: str, save_file: str,
             sent_max_len: int = 300):

        save_set = set()

        total_count = 0
        save_count = 0
        sent_over_count = 0
        with open(chain_file) as f, open(save_file, 'w') as wf:
            for line in f:
                l = line.strip().split('\t')

                if l == [ '' ]:
                    continue

                total_count += 1

                if total_count % 10000 == 0:
                    print(f'Processed: {total_count:,}, Saved: {save_count:,}.')

                pmid = l[ 0 ]
                sent_id = l[ 1 ]

                sent = l[ 2 ]
                chain = l[ 3: ]

                if (pmid, sent_id) in save_set:
                    continue
                save_set.add((pmid, sent_id))

                # 根据句子长度进行筛选
                if len(sent) > sent_max_len:
                    sent_over_count += 1
                    continue

                if not self.chain_filter(chain):
                    continue

                # todo: 根据基因，突变，调控词，生物过程在句子中的位置进行筛选
                save_count += 1
                wf.write(line)

        print(f'SentMaxLen: {sent_max_len:,}, FilterCount: {sent_over_count:,}.')
        print(f'SaveCount: {save_count:,} / {total_count:,}')
        print(f'{save_file} save done.')

    @staticmethod
    def save_no_norm_file(no_norm_count: dict, save_file: str):

        sorted_key = sorted(no_norm_count.keys(), key=lambda x: no_norm_count[ x ], reverse=True)

        with open(save_file, 'w') as wf:
            for _key in sorted_key:
                wf.write(f'{_key}\t{no_norm_count[ _key ]}\n')
        print(f'{save_file} save done.')

    def save_no_norm(self, save_path: str):

        no_norm_gene_save_file = f'{save_path}/no_norm_gene.tsv'
        no_norm_var_save_file = f'{save_path}/no_norm_var.tsv'
        no_norm_concept_save_file = f'{save_path}/no_norm_concept.tsv'

        self.save_no_norm_file(self.no_norm_gene, no_norm_gene_save_file)
        self.save_no_norm_file(self.no_norm_var, no_norm_var_save_file)
        self.save_no_norm_file(self.no_norm_concept, no_norm_concept_save_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-cf', dest='chain_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    parser.add_argument('-sl', dest='sent_len', type=int, default=300,
                        help='default: 300')
    args = parser.parse_args()

    chain_seeker = chain_selector()

    chain_seeker.main(args.chain_file, args.save_file, args.sent_len)

    base_path = '/'.join(os.path.abspath(args.save_file).split('/')[:-1])

    chain_seeker.save_no_norm(base_path)
