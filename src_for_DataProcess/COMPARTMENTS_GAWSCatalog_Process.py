# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 16/12/2022 15:32
@Author: XINZHI YAO
"""

import os
from collections import defaultdict

"""
改代码用于处理 COMPARTMENTS 和 GWASCatalog的数据 用于存到外链
"""

def read_compartments_file(input_file, symbol_set, int_file=False):

    matched_set = set()
    symbol_to_cc_set = defaultdict(set)
    with open(input_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            # print(l)
            # input()
            # not integrated file
            symbol = l[1]
            # print(symbol)
            cell_com = l[3]
            score = float(l[-1])

            if symbol not in symbol_set:
                continue
            matched_set.add(symbol)

            symbol_to_cc_set[symbol].add((cell_com, score))

    print(f'gene count: {len(symbol_to_cc_set.keys()):,}')
    print(f'Matched: {len(matched_set):,}/{len(symbol_set):,}')
    return symbol_to_cc_set

def read_symbol_file(symbol_file):
    symbol_set = set()

    with open(symbol_file) as f:
        for line in f:
            symbol_set.add(line.strip())
    print(f'Symbol: {len(symbol_set):,}')
    return symbol_set

def get_highest_score_compartment(symbol_to_cc_set: dict):

    symbol_to_highest_cc_set = defaultdict(set)
    for symbol, cc_set in symbol_to_cc_set.items():
        sort_cc_list = list(sorted(cc_set, key=lambda x: x[1], reverse=True))

        # (cc, score)
        highest_score = sort_cc_list[0][1]

        new_cc_set = set()
        for (compartment, score) in sort_cc_list:
            if score < highest_score:
                break
            symbol_to_highest_cc_set[symbol].add(compartment)
    return symbol_to_highest_cc_set

def save_symbol_to_cc_set(symbol_to_highest_cc_set: dict, save_file: str):

    max_len = 0
    with open(save_file, 'w') as wf:
        wf.write(f'Symbol\tCOMPARTMENT\n')
        for symbol, cc_set in symbol_to_highest_cc_set.items():
            cc_list = sorted(cc_set)
            cc_wf = ', '.join(cc_list)

            if len(cc_wf) > max_len:
                max_len = len(cc_wf)

            wf.write(f'{symbol}\t{cc_wf}\n')


    print(f'{save_file} saved, max_len: {max_len}.')


def main():

    symbol_total_file = '/mnt/disk1/xzyao/AD-PNRLE/data/external_data/symbol_total.tsv'

    symbol_set = read_symbol_file(symbol_total_file)

    compartments_file = '/mnt/disk1/xzyao/AD-PNRLE/data/external_data/COMPARTMENTS数据下载/human_compartment_knowledge_full.tsv'
    symbol_to_cc_set = read_compartments_file(compartments_file, symbol_set)

    symbol_to_highest_cc_set = get_highest_score_compartment(symbol_to_cc_set)

    save_file = '/mnt/disk1/xzyao/AD-PNRLE/data/external_data/symbol-compartment.tsv'

    save_symbol_to_cc_set(symbol_to_highest_cc_set, save_file)
    # compartments_file = '/mnt/disk1/xzyao/AD-PNRLE/data/external_data/COMPARTMENTS数据下载/human_compartment_integrated_full.tsv'
    # symbol_to_cc_set = read_compartments_file(compartments_file, symbol_set, True)






if __name__ == main():
    main()

