# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 22/09/2022 10:41
@Author: XINZHI YAO
"""

import os

import re

from tqdm import tqdm
import argparse
from collections import defaultdict


def read_mapping_file(mapping_file_path: str):
    total_mapping_file = f'{mapping_file_path}/go_hpo_mesh-umls.mapping.tsv'

    umls_to_concept = defaultdict(set)
    with open(total_mapping_file) as f:
        for line in f:
            cui, concept_id = line.strip().split('\t')

            umls_to_concept[cui].add(concept_id)

    return umls_to_concept

def id_extractor(concept: str):
    return re.findall(r"id: '(.*?)'", concept)[0]

def id_type(_id: str):
    if _id.startswith('GO:'):
        return 'GO'
    elif _id.startswith('HP:'):
        return 'HPO'
    elif _id.startswith('D'):
        return 'Mesh'
    else:
        return 'Unknown'

def main(iasis_path: str, save_file:str,
         mapping_file_path: str):

    umls_to_concept = read_mapping_file(mapping_file_path)

    iasis_file_list = os.listdir(iasis_path)

    wf = open(save_file, 'w')
    # fixme: new format of iASiS data: RelationType\tRelation\t
    wf.write(f'RelationType\tRelation\t'
             f'Source\tSourceID\tSourceType\tSourceSemanticType\t'
             f'Target\tTargetID\tTargetType\tTargetSemanticType\n')

    save_count = 0
    for iasis_file in tqdm(iasis_file_list):
        file_path = f'{iasis_path}/{iasis_file}'
        print(file_path)
        with open(file_path, encoding= u'utf-8',errors='ignore') as f:
            f.readline()
            for line in f:
                l = line.strip().split('\t')
                relation_type = l[0]
                relation = l[1]
                
                source_label = l[2]
                source_semantic_type = l[4]
                source = l[5]

                target_label = l[6]
                target_semantic_type = l[8]
                target = l[9]

                if source_label == 'None' or target_label == 'None':
                    continue

                source_id = id_extractor(source)
                source_id_set = umls_to_concept[source_id]

                target_id = id_extractor(target)
                target_id_set = umls_to_concept[target_id]

                for source_id in source_id_set:
                    for target_id in target_id_set:

                        source_type = id_type(source_id)
                        target_type = id_type(target_id)

                        save_count += 1
                        wf.write(f'{relation_type}\t{relation}\t'
                                 f'{source_label}\t{source_id}\t{source_type}\t'
                                 f'{source_semantic_type}\t'
                                 f'{target_label}\t{target_id}\t{target_type}\t'
                                 f'{target_semantic_type}\n')

    wf.close()
    print(f'{save_file} save done, {save_count:,} links saved.')


if __name__ == '__main__':
    parser =  argparse.ArgumentParser()
    parser.add_argument('-ip', dest='iasis_path', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)

    parser.add_argument('-mp', dest='mapping_file_path',
                        default='../result/go_hpo_mesh_umls_mapping',
                        help='default ../result/go_hpo_mesh_umls_mapping')
    args = parser.parse_args()

    main(args.iasis_path, args.save_file, args.mapping_file_path)


