# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 20/09/2022 20:21
@Author: XINZHI YAO
"""

import os
import argparse

from collections import defaultdict

def read_concept_file(concept_file: str):
    print(f'Reading {concept_file}.')
    concept_count = defaultdict(int)
    with open(concept_file) as f:
        f.readline()
        for line in f:
            concept, count = line.strip().split('\t')
            concept_count[concept] = int(count)
    return concept_count


def read_mapping_file(mapping_file: str):
    print(f'Reading {mapping_file}.')

    concept_to_umls = defaultdict(set)
    with open(mapping_file) as f:
        f.readline()
        for line in f:
            cui, concept = line.strip().split('\t')
            concept_to_umls[concept].add(cui)
    return concept_to_umls

def concept_convert(concept_count: dict, concept_umls_mapping: dict,
                    save_file: str):

    norm_count = 0
    with open(save_file, 'w') as wf:
        for concept, count in concept_count.items():
            if concept_umls_mapping.get(concept):
                norm_set = concept_umls_mapping[concept]
                norm_count += 1
                for norm in norm_set:
                    wf.write(f'{concept}\t{count}\t{norm}\n')
            else:
                wf.write(f'{concept}\t{count}\t-\n')
                continue
    print(f'{save_file} saved, {norm_count}/{len(concept_count.keys())} mapped.')


def main(concept_path: str, mapping_path: str, save_path: str):

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    # go concept/mapping/mapped file
    go_concept_file = f'{concept_path}/go.concept.tsv'
    go_umls_mapping_file = f'{mapping_path}/go_umls_mapping.tsv'
    go_mapped_file = f'{save_path}/go_umls_mapped.tsv'
    # hpo concept/mapping/mapped file
    hpo_concept_file = f'{concept_path}/hpo.concept.tsv'
    hpo_umls_mapping_file = f'{mapping_path}/hpo_umls_mapping.tsv'
    hpo_mapped_file = f'{save_path}/hpo_umls_mapped.tsv'
    # mesh concept/mapping/mapped file
    mesh_concept_file = f'{concept_path}/mesh.concept.tsv'
    mesh_umls_mapping_file = f'{mapping_path}/mesh_umls_mapping.tsv'
    mesh_mapped_file = f'{save_path}/mesh_umls_mapped.tsv'

    go_concept_count = read_concept_file(go_concept_file)
    hpo_concept_count = read_concept_file(hpo_concept_file)
    mesh_concept_count = read_concept_file(mesh_concept_file)

    go_to_umls = read_mapping_file(go_umls_mapping_file)
    hpo_to_umls = read_mapping_file(hpo_umls_mapping_file)
    mesh_to_umls = read_mapping_file(mesh_umls_mapping_file)

    print('GO-UMLS mapping.')
    concept_convert(go_concept_count, go_to_umls, go_mapped_file)
    print('HPO-UMLS mapping.')
    concept_convert(hpo_concept_count, hpo_to_umls, hpo_mapped_file)
    print('MeSH-UMLS mapping.')
    concept_convert(mesh_concept_count, mesh_to_umls, mesh_mapped_file)

    parser.add_argument('-cp', dest='concept_path', required=True)
    parser.add_argument('-mp', dest='mapping_path', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    args = parser.parse_args()

    main(args.concept_path, args.mapping_path, args.save_path)

