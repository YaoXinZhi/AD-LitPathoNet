# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 04/10/2022 20:16
@Author: XINZHI YAO
"""

import argparse
from collections import defaultdict
import os

def save_multi_id_file(concept_to_id: dict, save_file: str):

    with open(save_file, 'w') as wf:
        wf.write('Concept\tIDs\n')
        for concept, ids in concept_to_id.items():
            if len(ids) > 1:
                id_wf = '\t'.join(ids)
                wf.write(f'{concept}\t{id_wf}\n')
    print(f'{save_file} save done.')

def save_path_length_to_concept(path_length_to_concept: dict, save_file: str):

    with open(save_file, 'w') as wf:
        wf.write('PathLength\tEndConceptID\tEndConceptName\n')
        for path_length, concept_set in path_length_to_concept.items():
            for (concept_id, concept_name) in concept_set:
                wf.write(f'{path_length}\t{concept_id}\t{concept_name}\n')
    print(f'{save_file} save done.')

def save_concept_count_file(concept_count: dict, save_file: str, id_to_concept: dict):
    with open(save_file, 'w') as wf:
        wf.write('ID\tConcept\tCount\n')
        sort_id = sorted(concept_count.keys(),key = lambda x: concept_count[ x ], reverse = True)
        for _id in sort_id:
            wf.write(f'{_id}\t{id_to_concept[ _id ]}\t{concept_count[ _id ]}\n')
    print(f'{save_file} save done.')


def main(pathway_file: str, save_path: str, have_score: bool):
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    concept_count_file = f'{save_path}/concept_count.tsv'
    multi_id_concept_file = f'{save_path}/multi_id_concept.tsv'
    path_length_to_concept_file = f'{save_path}/concept_length_to_End-concept.tsv'

    path_len_to_concept = defaultdict(set)

    print(f'Processing {pathway_file}.')
    with open(pathway_file) as f:
        f.readline()
        concept_count = defaultdict(int)
        id_to_concept = {}
        while True:
            path_line = f.readline().strip()
            _ = f.readline().strip()
            label_line = f.readline().strip()
            line = f.readline()
            if not line:
                break

            pathway_split = path_line.split('\t')

            if have_score:
                path = [concept for concept in pathway_split[4:]]
                label_list = [label for label in label_line.split('\t')[4:]]
            else:
                path = [concept for concept in pathway_split[3:]]
                label_list = [label for label in label_line.split('\t')[3:]]
             
            
            
            if len(path) < 3 or len(label_list) < len(path):
                continue

            path_length = pathway_split[0]
            end_concept = path[-1]
            end_id = label_list[-1]

            path_len_to_concept[path_length].add((end_id, end_concept))

            for idx, concept in enumerate(path):
                label = label_list[idx]
                if label == '-':
                    continue
                concept_count[label] += 1
                id_to_concept[label] = concept

    concept_to_id = defaultdict(set)
    for _id, concept in id_to_concept.items():
        concept_to_id[ concept.lower() ].add(_id)

    save_concept_count_file(concept_count, concept_count_file, id_to_concept)
    save_multi_id_file(concept_to_id, multi_id_concept_file)
    save_path_length_to_concept(path_len_to_concept, path_length_to_concept_file)

if __name__ == '__main__':

    # statistic original pathway file (no score).
    parser = argparse.ArgumentParser()
    parser.add_argument('-pf', dest='pathway_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)

    parser.add_argument('-hs', dest='have_score',
                        default=False, action='store_true',
                        help='if have score in pathway file, default: False.')

    args = parser.parse_args()

    main(args.pathway_file, args.save_path, args.have_score)

