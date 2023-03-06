# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 13/10/2022 20:16
@Author: XINZHI YAO
"""

import argparse
import os


def main(pathway_path: str, save_path: str, chain_file: str, commend_file: str):
    # only pathway db
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    save_count = 0
    pathway_file_list = os.listdir(pathway_path)
    with open(commend_file, 'w') as wf:
        for pathway_file in pathway_file_list:

            if not (pathway_file.endswith('txt') or pathway_file.endswith('tsv')):
                continue

            prefix = '.'.join(pathway_file.split('.')[:2])

            file_path = f'{pathway_path}/{pathway_file}'

            commend = f'python3 iASiS_and_pathway_to_DB.py -pf {file_path} ' \
                      f'-cf {chain_file} -sp {save_path} -co pathway -pr {prefix}'

            wf.write(f'{commend}\n')
            save_count += 1

    print(f'{commend_file} save done, {save_count:,} commend line.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pp', dest='pathway_path', required=True)
    parser.add_argument('-cf', dest='chain_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    parser.add_argument('-cm', dest='commend_file',
                        default='bulk_pathway_to_db.sh',
                        help='default: bulk_pathway_to_db.sh')
    args = parser.parse_args()

    main(args.pathway_path, args.save_path, args.chain_file, args.commend_file)
