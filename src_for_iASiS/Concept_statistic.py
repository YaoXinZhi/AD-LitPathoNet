# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 20/09/2022 19:34
@Author: XINZHI YAO
"""

import os
import re

import argparse

from collections import defaultdict

def save_concept(save_file: str, concept_count: dict):

    sort_concept = sorted(concept_count, key=lambda x: concept_count[x], reverse=True)
    with open(save_file, 'w') as wf:
        wf.write(f'Concept\tCount\n')
        for concept in sort_concept:
            count = concept_count[concept]
            wf.write(f'{concept}\t{count}\n')
    print(f'{save_file} saved.')

def main(chain_file: str, save_path: str):
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    go_count = defaultdict(int)
    hpo_count = defaultdict(int)
    mesh_count = defaultdict(int)
    with open(chain_file) as f:
        for line in f:
            l = line.strip()

            go_list = re.findall(r'GO:\d+', l)
            hpo_list = re.findall(r'HP:\d+', l)
            mesh_list = re.findall(r'D\d+', l)

            for go in go_list:
                go_count[go] += 1
            for hpo in hpo_list:
                hpo_count[hpo] += 1
            for mesh in mesh_list:
                mesh_count[mesh] += 1

    go_save_file = f'{save_path}/go.concept.tsv'
    hpo_save_file = f'{save_path}/hpo.concept.tsv'
    mesh_save_file = f'{save_path}/mesh.concept.tsv'

    save_concept(go_save_file, go_count)
    save_concept(hpo_save_file, hpo_count)
    save_concept(mesh_save_file, mesh_count)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-cf', dest='chain_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    args = parser.parse_args()

    main(args.chain_file, args.save_path)
