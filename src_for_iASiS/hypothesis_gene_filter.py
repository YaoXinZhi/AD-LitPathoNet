# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 08/11/2022 21:18
@Author: XINZHI YAO
"""

import os
import argparse
from collections import defaultdict


"""
改代码用于根据假说基因来对pathway文件初步进行筛选

需要注意的是默认筛选的是没有filter的pathway文件
因此source列在第2列
筛选过后的在第4列

方便之后用linux命令进行筛选
"""

def main(pathway_file: str, save_path: str):

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    abeta_save_file = f'{save_path}/abeta.pathway.tsv'
    tau_save_file = f'{save_path}/tau.pathway.tsv'

    # abeta_hypo_gene = {'APP', 'BACE1', 'NCSTN', 'PSEN1', 'PSENEN', 'APH1A'}
    abeta_hypo_entrez = {'351', '23621', '23385', '5663', '55851', '51107'}
    tau_hypo_entrez = {'4137'}

    abeta_save_count = 0
    tau_save_count = 0
    process_count = 0
    with open(pathway_file) as f, open(abeta_save_file, 'w') as abeta_wf,\
        open(tau_save_file, 'w') as tau_wf:
        head = f.readline()
        head_split = head.strip().split('\t')
        source_col = head_split.index('Source')
        print(f'Source Column: {source_col}')

        abeta_wf.write(head)
        tau_wf.write(head)

        while True:
            path_line = f.readline().strip()
            type_line = f.readline().strip()
            label_line = f.readline().strip()

            process_count += 1
            if process_count % 50000 == 0:
                print(f'{process_count} pathway processed.')

            line = f.readline()
            if not line:
                break

            # print(path_line.split('\t')[source_col])
            source = path_line.split('\t')[source_col]

            if source in abeta_hypo_entrez:
                abeta_wf.write(path_line)
                abeta_wf.write(type_line)
                abeta_wf.write(label_line)
                abeta_wf.write('\n')
                abeta_save_count += 1
            if source in tau_hypo_entrez:
                tau_wf.write(path_line)
                tau_wf.write(type_line)
                tau_wf.write(label_line)
                tau_wf.write('\n')
                tau_save_count += 1

    print(f'{abeta_save_file} save done, {abeta_save_count:,} pathway saved.')
    print(f'{tau_save_file} save done, {tau_save_count:,} pathway saved.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pf', dest='pathway_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    args = parser.parse_args()

    main(args.pathway_file, args.save_path)
