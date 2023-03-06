# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 05/10/2022 10:08
@Author: XINZHI YAO
"""

import os
import re
import argparse


"""
该代码用于从PNRLE rich-evidence文件中找到现在pathway中PNRLE部分的证据
"""


def read_multi_id_file( multi_id_file: str):
    multi_id_mapping_dict = {}
    with open(multi_id_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            # concept = l[0]
            ids = l[ 1: ]

            selected_id = ''
            for _id in ids:
                if _id.startswith('GO'):
                    selected_id = _id
                    continue
                elif _id.startswith('HP') and not selected_id.startswith('GO'):
                    selected_id = _id
                elif (_id.startswith('D') or _id.startswith('C')) and not (
                        selected_id.startswith('GO') or selected_id.startswith('HP')):
                    selected_id = _id

            for _id in ids:
                multi_id_mapping_dict[ _id ] = selected_id

    return multi_id_mapping_dict

def read_evidence_file(evidence_db_file: str, multi_id_mapping_dict: dict):

    unique_edge_set = set()

    with open(evidence_db_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            unique_edge = eval(l[1])
            gene, var, rel, concept, theme_gene = unique_edge
            if multi_id_mapping_dict.get(concept):
                unique_edge_set.add((gene, var, rel, multi_id_mapping_dict[concept], theme_gene))
            else:
                unique_edge_set.add(unique_edge)

    print(f'{len(unique_edge_set)} unique edge have evidence.')
    return unique_edge_set


def main(pathway_file:str, evidence_db_file:str, multi_id_file: str):

    multi_id_mapping_dict = read_multi_id_file(multi_id_file)

    unique_edge_set = read_evidence_file(evidence_db_file, multi_id_mapping_dict)


    with open(pathway_file) as f:
        f.readline()
        while True:
            path_line = f.readline().strip()
            type_line = f.readline().strip()
            label_line = f.readline().strip()

            line = f.readline()
            if not line:
                break

            path_list = [elem for elem in path_line.split('\t')[5:] if elem != '-']

            label_list = [label for label in label_line.split('\t')[5:] if label != '-']

            relation = path_list[3]

            if relation not in {'Reg', 'PosReg', 'NegReg'}:
                print('not in relation')
                print(path_line)
                input()

            gene = label_list[0]
            var = label_list[1]
            concept = label_list[2]
            theme_gene = 'None'

            unique_edge = (gene, var, relation, concept, theme_gene)

            if unique_edge not in unique_edge_set:
                print('do not find evidence')
                print(unique_edge)
                print(path_line)
                input()





if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-ed', dest='evidence_db_file', required=True)

    parser.add_argument('-pf', dest='pathway_file', required=True)

    parser.add_argument('-mi', dest='multi_id_file', required=False,
                        default='../result/multi_id_concept.tsv',
                        help='default: ../result/multi_id_concept.tsv')
    args = parser.parse_args()

    main(args.pathway_file, args.evidence_db_file, args.multi_id_file)

