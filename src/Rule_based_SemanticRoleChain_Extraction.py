# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 25/04/2022 20:44
@Author: XINZHI YAO
"""

"""
该代码用于基于规则的语义角色链条抽取
1. 观测各个规则下生成的逻辑线条
2. 调整规则
"""

import re
import os
import argparse

from tqdm import tqdm

from collections import defaultdict


def add_theme(theme_to_entity: dict, main_chain_list: list):
    # add left ThemeOf
    left_theme = False
    continue_add = True
    chain_list = main_chain_list.copy()
    total_chain = chain_list.copy()
    while continue_add:
        continue_add = False

        for chain in chain_list:
            left_sub = chain[0]
            if theme_to_entity.get(left_sub):

                total_chain.remove(chain)
                for sub_info in theme_to_entity[ left_sub ]:
                    # Avoid ring
                    if sub_info not in chain:
                        continue_add = True
                        left_theme = True
                        total_chain.append([ sub_info, 'ThemeOf', *chain])
                    else:
                        continue_add = False

        chain_list = total_chain.copy()

    # add right ThemeOf
    right_theme = False
    continue_add = True
    while continue_add:
        continue_add = False
        for chain in chain_list:
            right_obj = chain[ -1 ]
            if theme_to_entity.get(right_obj):
                right_theme = True
                continue_add = True
                total_chain.remove(chain)

                for obj_info in theme_to_entity[right_obj]:
                    if obj_info not in chain:
                        continue_add = True
                        left_theme = True
                        total_chain.append([*chain, 'ThemeOf', obj_info])

        chain_list = total_chain.copy()

    return total_chain, left_theme, right_theme

def add_sub_events(sub_chain_list: list, theme_to_entity: dict):

    sub_chain_list_new = []
    ## delete first Reg for sub event chains
    for sub_chain in sub_chain_list:
        print(sub_chain)

        sub_reg, _, obj_reg = sub_chain
        for sub_entity in theme_to_entity[sub_reg]:
            sub_chain_list_new.append([sub_entity, 'CauseOf',obj_reg])

    sub_event_chains, left_theme, right_theme = add_theme(theme_to_entity, sub_chain_list_new)
    return sub_event_chains, left_theme, right_theme

def filter_both_themed(event_chains: str):

    reg_set = {'Reg', 'NegReg', 'PosReg'}

    # 筛选左右都添加了ThemeOf的链条
    # 判断为左边第一个节点和右边第一个节点都不是Reg
    filted_chain_list = []
    for chain in event_chains:
        left_entity_type = chain[0][1]
        right_entity_type = chain[-1][1]

        if left_entity_type not in reg_set \
                and right_entity_type not in reg_set:
            filted_chain_list.append(chain)
    return filted_chain_list

def rule_1(triple_set: set):
    reg_set = {'Reg', 'NegReg', 'PosReg'}
    # 1. find main/sub events chain
    # Var -- CauseOf -- Reg
    main_chain_list = []
    # Reg -- CauseOf -- Reg
    sub_chain_list = []
    for triple in triple_set:
        sub, sub_type, sub_offset, sub_norm,\
        obj, obj_type, obj_offset, obj_norm, relation = triple

        # fixme: detail of rule1
        if sub_type in {'Var'} \
            and obj_type in reg_set \
            and relation == 'CauseOf':

            main_chain_list.append([(sub, sub_type, sub_offset, sub_norm), 'CauseOf', (obj, obj_type, obj_offset, obj_norm)])

        # todo: Sub event Reg -- CauseOf Reg
        if sub_type in reg_set \
            and obj_type in reg_set \
            and relation == 'CauseOf':
            sub_chain_list.append([(sub, sub_type, sub_offset, sub_norm), 'CauseOf', (obj, obj_type, obj_offset, obj_norm)])

    # 2. ThemeOf Statistics
    theme_to_entity = defaultdict(set)
    for triple in triple_set:
        sub, sub_type, sub_offset, sub_norm,\
        obj, obj_type, obj_offset, obj_norm, relation = triple

        # fixme: detail of theme
        if relation == 'ThemeOf':
            # entity_to_theme[(sub, sub_type, sub_offset)].add((obj, obj_type, obj_offset))
            theme_to_entity[(obj, obj_type, obj_offset, obj_norm)].add((sub, sub_type, sub_offset, sub_norm))

    # 3. Add ThemeOf for main events chain
    main_event_chain, _, _ = add_theme(theme_to_entity, main_chain_list)

    # 4. Add ThemeOf for sub event chains
    ## 4.1 version 1. delete first reg
    # if sub_chain_list:
    #     sub_event_chain, _, _ = add_sub_events(sub_chain_list, theme_to_entity)
    ## 4.2 retain the first reg
    sub_event_chain, _, _ = add_theme(theme_to_entity, sub_chain_list)

    return main_event_chain, sub_event_chain, _, _


def read_infer_file(infer_file: str):

    """
    Entity1 Entity1 Type    Entity1 Position        Entity1 Tool-Tag        Entity1 AGAC    Entity1 NormEntity1 Norm-Count
    Entity2 Entity2 Type    Entity2 Position        Entity2 Tool-Tag        Entity2AGAC    Entity2 Norm    Entity2 Norm-Count
    None
    Sentence        PMID    Sentence ID     Relation
    """

    print(f'Reading {infer_file}.')
    sent_id_to_infer_set = defaultdict(set)

    with open(infer_file) as f:
        head_line = f.readline()
        for line in f:
            l = line.strip().split('\t')
            relation = l[ -1 ]
            sent_id = l[-2]
            pmid = l[-3]
            sent = l[-4]
            if relation != 'NoRelation':
                sent_id_to_infer_set[(pmid, sent_id, sent)].add(line)

    return sent_id_to_infer_set, head_line

def infer_set_to_triple_set(infer_set: set):

    triple_set = set()
    for line in infer_set:
        l = line.strip().split('\t')

        entity1 = l[ 0 ]
        entity1_type = l[ 1 ]
        entity1_offset = l[ 2 ]
        entity1_norm = l[5]

        entity2 = l[ 7 ]
        entity2_type = l[ 8 ]
        entity2_offset = l[ 9 ]
        entity2_norm = l[12]

        # sent = l[-4]
        # pmid = l[-3]
        # sent_id = l[ -2 ]
        relation = l[-1]

        if relation != 'NoRelation' \
                and (entity1, entity1_type, entity1_offset) != (entity2, entity2_type, entity2_offset):
            triple_set.add((entity1, entity1_type, entity1_offset, entity1_norm,
                            entity2, entity2_type, entity2_offset, entity2_norm,
                            relation))
    return triple_set


def chain_filter(event_chains: list, no_annotation_go_set: set,
                 entrez_to_rs: dict, rs_to_entrez: dict,
                 rs_entrez_to_symbol: dict):
    """
    rule 1:
    1. 2022-05-26 ThemeOf空缺的GO，如果comment不建议标注，则去除
    rule 2:
    2. 2022-05-26 如果chain开头为Gene ThemeOf rs-Var, 且该rs_id 和entrez没有对应上 则去除
    3. 2022-06-06 如果开头为Gene ThemeOf rs-Var, 如果rs_id和entrez没有对应上 则重新给他加上对应的个呢
    rule 3:
    4. 2022-06-06 如果开头为 rs-Var, 则加上对应的Gene
    """
    saved_chain_list = []
    for chain in event_chains:
        # filter rule 1.
        left_element = chain[-1]
        element_id = left_element[-1].split('-')[-1]
        if element_id.startswith('GO') and element_id in no_annotation_go_set:
            continue

        # filter rule 2
        if chain[1] == 'ThemeOf' and entrez_to_rs:
            first = chain[0]
            second = chain[2]
            # var 以rs开头
            if (first[1] == 'Gene' and first[-1] != '-') \
                    and (second[1] == 'Var' and second[-1] != '-')\
                    and re.findall(r'[Rr]s[0-9]+', second[-1]):
                entrez = first[-1].split('-')[-1]
                rs_id = second[-1].split('-')[-1]

                if rs_id not in entrez_to_rs[entrez] \
                    and not rs_to_entrez.get(rs_id):
                    continue
                cor_entrez = rs_to_entrez[rs_id]
                symbol = rs_entrez_to_symbol[cor_entrez]
                cor_gene = (symbol, 'Gene', 'dbSNP', f'{symbol}-Gene-{cor_entrez}')
                chain[0] = cor_gene


        # filter rule 3
        if chain[0][1] == 'Var' and chain[0][-1] != '-':
            # rs-Var is first element in chain
            if re.findall(r'[Rr]s[0-9]+', chain[0][1]):

                rs_id = chain[0][-1].split('-')[-1]

                if not rs_to_entrez.get(rs_id):
                    continue
                entrez = rs_to_entrez[rs_id]
                symbol = rs_entrez_to_symbol[entrez]
                gene = (symbol, 'Gene', 'dbSNP', f'{symbol}-Gene-{entrez}')
                chain.insert(0, gene)
            # The first element is No-normalized variation.
            else:
                continue

        saved_chain_list.append(chain)

    return saved_chain_list

def read_obo(obo_file: str):
    id_to_name = {}

    # for MF, BP, CC
    id_to_type = {}

    _id = ''

    alt_id_set = set()
    no_annotation_id_set = set()
    with open(obo_file) as f:
        for line in f:

            l = line.strip().split(': ')
            if len(l) < 2:
                continue
            if line.startswith('id:'):
                _id = l[ 1 ]
                alt_id_set = set()
            if line.startswith('name:'):
                name = l[ 1 ]
                id_to_name[ _id ] = name
            if line.startswith('comment:'):
                if 'not be used for direct annotation' in line:
                    no_annotation_id_set.add(_id)
                    no_annotation_id_set.update(alt_id_set)
            if line.startswith('subset:'):
                if 'gocheck_do_not_annotate' in line or 'gocheck_do_not_manually_annotate' in line:
                    no_annotation_id_set.add(_id)
                    no_annotation_id_set.update(alt_id_set)
            if line.startswith('namespace:'):
                term_type = l[ 1 ]
                id_to_type[ _id ] = term_type
            if line.startswith('alt_id:'):
                alt_id = l[ 1 ]
                alt_id_set.add(alt_id)
                id_to_name[ alt_id ] = name

    return id_to_name, no_annotation_id_set, id_to_type

def read_rs_file(rs_file: str):

    rs_to_entrez = {}
    entrez_to_rs = defaultdict(set)
    rs_entrez_to_symbol = {}
    with open(rs_file) as f:
        f.readline()
        for line in f:
            rs_id, entrez, symbol = line.strip().split('\t')

            entrez_to_rs[entrez].add(rs_id)
            rs_to_entrez[rs_id] = entrez
            rs_entrez_to_symbol[entrez] = symbol

    return entrez_to_rs, rs_to_entrez, rs_entrez_to_symbol


def complete_infer_chain_extraction(infer_file: str, save_file: str,
                                    only_themed: bool, only_chain: bool,
                                    go_file: str, rs_file: str):

    sent_id_to_infer_set, head_line = read_infer_file(infer_file)


    _, no_annotation_go_set, _ = read_obo(go_file)


    if rs_file is None:
        entrez_to_rs = defaultdict(set)
        rs_to_entrez = {}
        rs_entrez_to_symbol = {}
    else:
        entrez_to_rs, rs_to_entrez, rs_entrez_to_symbol = read_rs_file(rs_file)

    wf = open(save_file, 'w')
    if not only_chain:
        wf.write(head_line)
    for (pmid, sent_id, sent), infer_set in sent_id_to_infer_set.items():
        triple_set = infer_set_to_triple_set(infer_set)

        if len(infer_set) < 2:
            continue

        # rule 1
        main_event_chains, sub_event_chains, left_theme, right_theme = rule_1(triple_set)

        if only_themed:
            main_event_chains = filter_both_themed(main_event_chains)
            sub_event_chains = filter_both_themed(sub_event_chains)

        main_event_chains = chain_filter(main_event_chains, no_annotation_go_set, entrez_to_rs, rs_to_entrez, rs_entrez_to_symbol)
        sub_event_chains = chain_filter(sub_event_chains, no_annotation_go_set, entrez_to_rs, rs_to_entrez, rs_entrez_to_symbol)

        if main_event_chains or sub_event_chains:
            # save result
            if not only_chain:
                for infer_line in infer_set:
                    wf.write(infer_line)
                for main_idx, chain in enumerate(main_event_chains):
                    chain_wf = '\t'.join(map(str, chain))
                    wf.write(f'Main-SRC-{main_idx}:\t{chain_wf}\n')
                for sub_idx, chain in enumerate(sub_event_chains):
                    chain_wf = '\t'.join(map(str, chain))
                    wf.write(f'Sub-SRC-{sub_idx}:\t{chain_wf}\n')
            else:
                for main_idx, chain in enumerate(main_event_chains):
                    chain_wf = '\t'.join(map(str, chain))
                    wf.write(f'{pmid}\t{sent_id}\t{sent}\t{chain_wf}\n')
                for sub_idx, chain in enumerate(sub_event_chains):
                    chain_wf = '\t'.join(map(str, chain))
                    wf.write(f'{pmid}\t{sent_id}\t{sent}\t{chain_wf}\n')
            wf.write('\n')
            # wf.write(f'Left theme: {left_theme}\tRight theme:{right_theme}\n\n')
    wf.close()
    print(f'{infer_file} processed, {save_file} saved.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rule based Semantic Role Chains Extraction.')
    parser.add_argument('--infer_file', dest='infer_file',
                        help='infer file with entity_norm.')
    parser.add_argument('--save_file', dest='save_file',
                        help='Save file for Extracted Semantic Role Chains')

    parser.add_argument('--only_themed', dest='only_themed', action='store_true',
                        default=False,
                        help='Only the chains with ThemeOf added to both left and right are retained.')

    parser.add_argument('--only_chain', dest='only_chain',
                        action='store_true', default=False,
                        help='save only chain in save_file, default: False.')

    parser.add_argument('--read_all', dest='read_all',
                        action='store_true', default=False,
                        help='Read in all infer files to avoid wrong sent_id.')

    parser.add_argument('--go_file', dest='go_file',
                        default='../data/obo-data/go.obo',
                        help='default: ../data/obo-data/go.obo.')

    # entrez-rsID file for rs mutations filter.
    parser.add_argument('--rs_file', dest='rs_file',
                        default=None,
                        help='rs entrez file for ThemeOf filter.')

    args = parser.parse_args()

    print('-'*60)
    if args.read_all:
        complete_infer_chain_extraction(args.infer_file, args.save_file,
                                        args.only_themed, args.only_chain,
                                        args.go_file, args.rs_file)
    # else:
    #     chain_extraction(args.infer_file, args.save_file, args.only_themed,
    #                      args.only_chain, args.go_file)
    print('-'*60)


