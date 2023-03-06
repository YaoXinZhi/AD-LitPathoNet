# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 31/07/2022 10:00
@Author: yao
"""

import os

import argparse
from collections import defaultdict


# 判断两个offset重叠
def is_offset_overleap(offset1: tuple, offset2: tuple):
    start1, end1 = map(int, offset1)
    start2, end2 = map(int, offset2)

    if start1 <= start2 <= end1 or start1 <= end2 <= end1:
        return True
    else:
        return False


def read_tagging_file(tagging_file: str):
    sent_id_to_tag = defaultdict(set)
    sent_id_to_agac = defaultdict(set)
    with open(tagging_file) as f:
        for line in f:
            l = line.strip().split('\t')

            sent_id = l[ 0 ]

            tag_list = [ eval(tag) for tag in l[ 3: ] ]

            for tag in tag_list:
                if len(tag) == 4 and tag[ 3 ] != 'None':
                    sent_id_to_tag[ sent_id ].add(tag)
                if len(tag) == 3:
                    sent_id_to_agac[ sent_id ].add(tag)
    return sent_id_to_tag, sent_id_to_agac


# 用于合并infer文件和tagging文件
def tagging_infer_merge(infer_file: str, sent_id_to_tag: dict, sent_id_to_agac: dict, save_file):
    wf = open(save_file, 'w')

    wf.write(
        f'Entity1\tEntity1 Type\tEntity1 Position\tEntity1 Tool-Tag\tEntity1 AGAC\tEntity1 Norm\tEntity1 Norm-Count\tEntity2\tEntity2 Type\tEntity2 Position\tEntity2 Tool-Tag\tEntity2 AGAC\tEntity2 Norm\tEntity2 Norm-Count\tNone\tSentence\tPMID\tSentence ID\tRelation\n')
    with open(infer_file) as f:
        for line in f:
            l = line.strip().split('\t')

            sent_id = l[ -2 ]

            entity1 = l[ 0 ]
            entity1_offset = eval(l[ 2 ])
            entity2 = l[ 3 ]
            entity2_offset = eval(l[ 5 ])

            # normalized tags for write
            tag_set = sent_id_to_tag[ sent_id ]

            entity1_norm = set()
            entity2_norm = set()
            entity2_candidate_norm = set()
            entity1_candidate_norm = set()
            for tag in tag_set:
                # 添加潜在标准化节点
                mention, _type, entity_id, position = tag

                # 这是看其他工具对这个有重叠的实体有没有标注 方便实体标注化
                entity1_is_overleap = is_offset_overleap(entity1_offset, position)
                if entity1_is_overleap:
                    if entity_id not in {'None', '-'}:
                        entity1_candidate_norm.add(f'{mention}-{_type}-{entity_id}')

                entity2_is_overleap = is_offset_overleap(entity2_offset, position)
                if entity2_is_overleap:
                    if entity_id not in {'None', '-'}:
                        entity2_candidate_norm.add(f'{mention}-{_type}-{entity_id}')

                # 这是判断其他工具对这个完全一样的mention有没有其他标注
                if mention == entity1:
                    entity1_norm.add(f'{_type}-{entity_id}')
                    if entity_id not in {'None', '-'}:
                        entity1_candidate_norm.add(f'{mention}-{_type}-{entity_id}')
                if mention == entity2:
                    entity2_norm.add(f'{_type}-{entity_id}')
                    if entity_id not in {'None', '-'}:
                        entity2_candidate_norm.add(f'{mention}-{_type}-{entity_id}')

            if entity1_norm:
                entity1_norm_wf = ', '.join(entity1_norm)
            else:
                entity1_norm_wf = '-'

            if entity2_norm:
                entity2_norm_wf = ', '.join(entity2_norm)
            else:
                entity2_norm_wf = '-'

            l.insert(6, entity2_norm_wf)
            l.insert(3, entity1_norm_wf)

            if entity1_candidate_norm:
                entity1_candidate_norm_wf = ', '.join(entity1_candidate_norm)
            else:
                entity1_candidate_norm_wf = '-'

            if entity2_candidate_norm:
                entity2_candidate_norm_wf = ', '.join(entity2_candidate_norm)
            else:
                entity2_candidate_norm_wf = '-'
            # print(entity2_norm)
            # print(entity2_candidate_norm)
            # print(len(entity2_candidate_norm))
            # print(entity2_candidate_norm_wf)
            # print(entity1_candidate_norm_wf)
            # print(l)
            # input()
            l.insert(8, entity2_candidate_norm_wf)
            # print(l)
            # input()
            l.insert(9, str(len(entity2_candidate_norm)))
            # print(l)
            # input()
            l.insert(4, entity1_candidate_norm_wf)
            # print(l)
            # input()
            l.insert(5, str(len(entity1_candidate_norm)))
            # print(l)
            # input()

            # AGAC tags
            agac_tag_set = sent_id_to_agac[ sent_id ]

            entity1_agac = set()
            entity2_agac = set()
            for tag in agac_tag_set:
                mention, _type, _ = tag
                if mention == entity1:
                    # print(mention, _type)
                    entity1_agac.add(_type)
                if mention == entity2:
                    # print(mention, _type)
                    entity2_agac.add(_type)

            if entity1_agac:
                entity1_agac_wf = ', '.join(entity1_agac)
            else:
                entity1_agac_wf = '-'

            if entity2_agac:
                entity2_agac_wf = ', '.join(entity2_agac)
            else:
                entity2_agac_wf = '-'
            # print(l)
            # input()
            l.insert(10, entity2_agac_wf)
            # print(l)
            # input()
            l.insert(4, entity1_agac_wf)
            # print(l)
            # input()

            line_wf = '\t'.join(l)

            wf.write(f'{line_wf}\n')

    wf.close()
    print(f'{save_file} save done.')

def main(tagging_file: str, infer_file: str, save_file: str):

    sent_id_to_tag, sent_id_to_agac = read_tagging_file(tagging_file)

    tagging_infer_merge(infer_file, sent_id_to_tag, sent_id_to_agac, save_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tagging Infer Merge.')
    parser.add_argument('-tg', dest='tagging_file', required=True)
    parser.add_argument('-if', dest='infer_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    main(args.tagging_file, args.infer_file, args.save_file)



