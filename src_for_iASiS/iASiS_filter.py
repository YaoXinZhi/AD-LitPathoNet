# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 26/09/2022 16:41
@Author: XINZHI YAO
"""

import os

import argparse

def main(input_file: str, save_file: str):
    functionally_relations = ['AFFECTS', 'MANAGES', 'TREATS', 'DISRUPTS', 'COMPLICATES', 'INTERACTS_WITH', 'PREVENTS',
             'BRINGS_ABOUT', 'PRODUCES', 'CAUSES', 'PERFORMS', 'CARRIES_OUT', 'EXHIBITS',
             'PRACTICES', 'OCCURS_IN', 'PROCESS_OF', 'USES', 'MANIFESTATION_OF', 'INDICATES',
             'RESULT_OF']

    with open(input_file) as f, open(save_file, 'w') as wf:
        for line in f:
            l = line.strip().split('\t')
            semantic_relation = l[0]

            save_bool = False
            for rel in functionally_relations:
                if semantic_relation.startswith(rel):
                    save_bool = True

            if save_bool:
                wf.write(line)
    print(f'{save_file} saved.')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-if', dest='input_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()


    main(args.input_file, args.save_file)

