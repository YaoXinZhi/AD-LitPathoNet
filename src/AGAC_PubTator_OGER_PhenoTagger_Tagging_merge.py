# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 03/08/2022 15:07
@Author: yao
"""
import argparse
import os

from collections import defaultdict
from itertools import permutations

# read HPO file
def read_hpo_tagger_file(hpo_tag_file: str):

    print('Reading HPO tagger file.')
    pmid_to_hpo = defaultdict(set)
    with open(hpo_tag_file) as f:
        for line in f:
            l = line.strip().split('\t')
            if len(l) < 5:
                continue

            pmid = l[0]

            tag = l[3]
            tag_type = l[4]
            tag_id = l[5]
            # confidence = l[6]

            pmid_to_hpo[pmid].add((tag, tag_type, tag_id))

    return pmid_to_hpo

# read OGER file
def read_oger_tagger_file(oger_tag_file: str):
    pmid_sent_id_to_tag = {}

    with open(oger_tag_file) as f:
        for line in f:
            l = line.strip().split('\t')

            sent_id = l[ 0 ]
            pmid = l[ 1 ]
            sent = l[ 2 ]
            tag = eval(l[ 3 ])

            pmid_sent_id_to_tag[ (pmid, sent_id) ] = tag
    return pmid_sent_id_to_tag

def find_offset(token: str, sent: str):
    offset_set = set()
    token_end = 0
    while True:
        token_start = sent.find(token, token_end)
        token_end = token_start + len(token)
        if token_start == -1:
            break
        offset_set.add((token_start, token_end))

    return offset_set

# 将AGAC，PubTator，OGER 和 PhenoTagger的标签合并成一个tagging文件
def merge_tag(tagging_file: str, pmid_to_hpo: dict,
              pmid_sent_id_to_tag: dict, save_file: str):

    wf = open(save_file, 'w')
    with open(tagging_file) as f:
        for line in f:
            l = line.strip().split('\t')

            # PMC:0001
            sent_idx = l[ 0 ]
            pmid = l[ 1 ]
            sentence = l[ 2 ]

            tag_set = set()
            for tag in l[ 3: ]:
                tag_set.add(eval(tag))

            # add HPO tag
            if pmid_to_hpo.get(pmid):
                for hpo_tag in pmid_to_hpo[ pmid ]:
                    tag = hpo_tag[ 0 ]
                    if tag in sentence:
                        mention, hpo_type, hpo_id = hpo_tag

                        offset_set = find_offset(mention, sentence)

                        for offset in offset_set:
                            tag_set.add((*hpo_tag, offset))
            # add GO tag
            sent_id = sent_idx.split(':')[ 1 ]
            if pmid_sent_id_to_tag.get((pmid, sent_id)):

                for go_tag in pmid_sent_id_to_tag[ (pmid, sent_id) ]:
                    mention = go_tag[ 0 ]
                    go_type = go_tag[ 2 ]
                    go_id = go_tag[ 3 ]
                    offset = go_tag[ 4 ]

                    tag_set.add((mention, go_type, go_id, offset))

            tag_wf = '\t'.join(map(str, tag_set))
            wf.write(f'{sent_idx}\t{pmid}\t{sentence}\t{tag_wf}\n')
    wf.close()
    print(f'{save_file} save done.')


# tagging merge
def tagging_merge(pubtator_set: set, model_tag_set: set):
    # Disease, Interaction, Enzyme, NegReg, Reg, Protein, MPA, CPA, Var,
    # todo: 20220414 把OGER和PhenoTagger的标签也加进去
    pubtator_label_to_agac = {
        # Pubtator
        'Disease': 'Disease',
        'DNAMutation': 'Var',
        'Gene': 'Gene',
        'ProteinMutation': 'Var',
        'SNP': 'Var',
        # SETH
        'SUBSTITUTION': 'Var',
        'DBSNP_MENTION': 'Var',
        'DELETION': 'Var',
        'FRAMESHIFT': 'Var',
        'INSERTION': 'Var',
        # PhenoTagger
        'Phenotype': 'Disease',
        # OGER
        # 这一类别的GO不太好归进去
        # 'cellular_component': '',
        'molecular_function': 'MPA',
        'biological_process': 'CPA',
    }

    special_save = {'p.M239V', 'AD'}

    pubtator_tag_set = {(tag[ 0 ], pubtator_label_to_agac[ tag[ 1 ] ], tag[ -1 ]) for tag in pubtator_set
                        if pubtator_label_to_agac.get(tag[ 1 ])}

    pubtator_token2label = {token: label for token, label, _ in pubtator_tag_set}

    # 2021 05 25
    # if a token appears in both pubtator and model labels,
    # but has a different label,
    # keep pubtator label.
    # Avoid nested entities
    merge_tag_set = pubtator_tag_set.copy()

    for m_token, m_label, m_offset in model_tag_set:

        (m_start, m_end) = m_offset

        m_start = int(m_start)
        m_end = int(m_end)

        save_flag = True
        for p_token, p_label, p_offset in pubtator_tag_set:
            (p_start, p_end) = p_offset

            p_start = int(p_start)
            p_end = int(p_end)
            # m_token in p_token
            if m_start >= p_start and m_end <= p_end:
                save_flag = False
                continue
            # p_token in m_token
            if p_start >= m_start and p_end <= m_end:
                # do not replace 527
                # if (p_token, p_label, p_offset) in merge_tag_set:
                #     merge_tag_set.remove((p_token, p_label, p_offset))
                # merge_tag_set.add((m_token, m_label, m_offset))
                save_flag = False
                continue

            # m_token cross p_token and m_token is in front of p_token
            if m_start <= p_start and (p_start <= m_end <= m_end):
                save_flag = False
                continue

            # p_token cross m_token and p_token is in front of m_token
            if p_start <= m_start and (m_start <= p_end <= m_end):
                save_flag = False
                pass

        if save_flag:
            merge_tag_set.add((m_token, m_label, m_offset))

    # merge label to pubtator label
    merge_tag_set_cp = merge_tag_set.copy()
    for token, label, offset in merge_tag_set_cp:
        if pubtator_token2label.get(token):
            if label == pubtator_token2label[ token ]:
                continue
            else:
                merge_tag_set.remove((token, label, offset))
                merge_tag_set.add((token, pubtator_token2label[ token ], offset))

    merge_tag_set_cp = merge_tag_set.copy()

    # tagging filter
    for info in merge_tag_set_cp:
        (token, label, (start, end)) = info
        if token in special_save:
            continue
        if len(token) <= 2:
            if info in merge_tag_set:
                merge_tag_set.remove(info)
                continue

        if '(' == token[ 0 ] \
                and ')' == token[ -1 ] \
                and label in [ 'Protein', 'Gene', 'Var' ]:
            merge_tag_set.remove(info)
            merge_tag_set.add((token[ 1:-1 ], label, (start + 1, end - 1)))
            continue

        if '(' == token[ 0 ] \
                and label in [ 'Protein', 'Gene', 'Var' ] \
                and len(token) >= 3:
            merge_tag_set.remove(info)
            merge_tag_set.add((token[ 1: ], label, (start + 1, end)))
            continue

        if ')' == token[ -1 ] and label in [ 'Protein', 'Gene', 'Var' ] and len(token) >= 3:
            merge_tag_set.remove(info)
            merge_tag_set.add((token[ :-1 ], label, (start, end - 1)))
            continue

        if ';' == token[ -1 ] or ':' == token[ -1 ]:
            merge_tag_set.remove(info)
            token = token[ :-1 ]
            end -= 1
            info = (token, label, (start, end))
            merge_tag_set.add(info)
            continue

        if "'s" == token[ -2: ]:
            merge_tag_set.remove(info)
            token = token[ :-2 ]
            end -= 2
            info = (token, label, (start, end))
            merge_tag_set.add(info)
            continue

        if '(' in token \
                or ')' in token \
                or ',' in token \
                or ' and ' in token \
                or len(token) > 50:
            if '(' in token \
                    and ')' in token \
                    and label in [ 'Protein', 'Gene', 'Var' ] \
                    and 'and' not in token \
                    and 3 <= len(token) <= 30:
                continue
            merge_tag_set.remove(info)
            continue

    merge_tag_set_cp = merge_tag_set.copy()
    # tagging merge
    for info_1 in merge_tag_set_cp:
        for info_2 in merge_tag_set_cp:
            (token_1, label_1, (start_1, end_1)) = info_1
            (token_2, label_2, (start_2, end_2)) = info_2

            start_1 = int(start_1)
            end_1 = int(end_1)
            start_2 = int(start_2)
            end_2 = int(end_2)

            if token_1 in special_save:
                continue

            if len(token_1.split()) == 1 and label_1 in [ 'Protein', 'Gene', 'Var', 'Disease' ]:
                continue

            if (start_1 >= start_2 and end_1 <= end_2) \
                    and info_1 != info_2:

                if start_1 == start_2 and end_1 == end_2:
                    if label_1 == 'Protein' and label_2 == 'Gene':
                        if info_1 in merge_tag_set:
                            merge_tag_set.remove(info_1)
                    if label_2 == 'Protein' and label_1 == 'Gene':
                        if info_2 in merge_tag_set:
                            merge_tag_set.remove(info_2)
                else:
                    if info_1 in merge_tag_set:
                        # print(info_1)
                        merge_tag_set.remove(info_1)
                        # merge_tag_set.remove((token_1, label_1, (start_1, end_1)))
    return merge_tag_set


# 用tagging文件生成 任务二的输入文件
# 生成的规则和原来 NER infer_result_process.py 里用的一样
def generate_re_input_file_from_tagging(tagging_file: str, save_file: str):
    non_self_label = {'Protein', 'Gene', "Enzyme", 'Var', 'Disease', 'MPA', 'CPA'}

    reg_set = {'Reg', 'PosReg', 'NegReg'}

    save_count = 0
    wf = open(save_file, 'w')
    with open(tagging_file) as f:
        for line in f:
            l = line.strip().split('\t')
            sent_id = l[ 0 ]
            pmid = l[ 1 ]
            sentence = l[ 2 ]
            # print(l)
            tag_set = [ eval(tag) for tag in l[ 3: ] ]

            # print(tag_set)
            # print(len(tag_set))
            # input()
            agac_tag_set = set()
            other_tag_set = set()

            # 分AGAC标签和标准化了的（PubTator，OGER，PhenoTagger）的标签
            for tag in tag_set:
                if len(tag) == 3:
                    agac_tag_set.add(tag)
                else:
                    other_tag_set.add(tag)
            # print(agac_tag_set)
            # print(f'agac: {len(agac_tag_set)}')
            # print(other_tag_set)
            # print(f'other: {len(other_tag_set)}')
            # input()

            merged_tagging_set = tagging_merge(other_tag_set, agac_tag_set)

            # print(merged_tagging_set)
            # print(len(merged_tagging_set))
            # input()

            for (tag1, tag2) in permutations(merged_tagging_set, 2):

                token1, label1, offset1 = tag1
                token2, label2, offset2 = tag2

                # if  token1 == 'growth' or token2 == 'growth':
                # if token2 == 'growth':
                #     print(token1, token2)
                #     input()

                if label1 in non_self_label and label2 in non_self_label and label1 == label2:
                    continue
                if label1 in reg_set and label2 not in reg_set:
                    # print(label1, label2)
                    # input('continue')
                    continue
                save_count += 1

                # if token2 == 'growth':
                #     print('writed')
                #     input()
                wf.write(f'{token1}\t{label1}\t{offset1}\t'
                         f'{token2}\t{label2}\t{offset2}\t'
                         f'None\t{sentence}\t{pmid}\t{sent_id}\n')

    wf.close()
    print(f'{save_file} save done.')

def main(tagging_file: str, oger_file: str, phenotagger_file: str,
         merge_save_file: str, re_input_save_file: str):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='AGAC PubTator OGER PhenoTagger Tagging Merger.')
    parser.add_argument('--action', dest='action', choices=['merge_tag', 'generate_re_input', 'both'],
                        default='both')

    parser.add_argument('-tf', dest='tagging_file',
                        required=False, default=None)
    parser.add_argument('-of', dest='oger_file',
                        required=False, default=None)
    parser.add_argument('-pf', dest='phenotagger_file',
                        required=False, default=None)
    parser.add_argument('-ts', dest='tagging_save_file',
                        required=False, default=None)
    parser.add_argument('-rs', dest='re_input_save_file',
                        required=False, default=None)
    args = parser.parse_args()

    if args.action in ['merge_tag', 'both']:
        
        if args.tagging_file is None \
            or args.oger_file is None \
            or args.phenotagger_file is None \
            or args.tagging_save_file is None:
            raise ValueError('parameters is None.')
        
        pmid_to_hpo = read_hpo_tagger_file(args.phenotagger_file)

        pmid_sent_id_to_tag = read_oger_tagger_file(args.oger_file)

        merge_tag(args.tagging_file, pmid_to_hpo, pmid_sent_id_to_tag, args.tagging_save_file)

    if args.action in ['generate_re_input', 'both']:

        if args.re_input_save_file is None \
            or args.tagging_save_file is None:
            raise ValueError('parameters is None.')

        generate_re_input_file_from_tagging(args.tagging_save_file, args.re_input_save_file)






