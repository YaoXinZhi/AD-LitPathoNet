# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 10/10/2022 16:00
@Author: XINZHI YAO
"""

"""
改代码用于从全pathway文件中拿出通过linux命令行找出的假说通路
对应成三行的pathway格式
"""

import os
import argparse

def main(pathway_file: str, hypothesis_path: str, match_prefix: str, save_file: str):

    file_list = [_file for _file in os.listdir(hypothesis_path) if _file.startswith(match_prefix)]

    # read selected pathway
    selected_pathway_set = set()
    for _file in file_list:
        file_path = f'{hypothesis_path}/{_file}'
        with open(file_path) as f:
            for line in f:
                l = line.strip()
                selected_pathway_set.add(l)
    print(f'{len(selected_pathway_set):,} pathway selected.')

    # pathway pick up
    save_count = 0
    with open(pathway_file) as f, open(save_file, 'w') as wf:
        f.readline()
        wf.write('ChainIdx\tScore\tPathLength\tSource\tTarget\tChain\n')
        while True:
            path_line = f.readline().strip()
            type_line = f.readline().strip()
            label_line = f.readline().strip()
            line = f.readline()
            if not line:
                break

            if path_line in selected_pathway_set:
                save_count += 1
                wf.write(f'{path_line}\n{type_line}\n{label_line}\n')
                wf.write('\n')

    print(f'{save_file} save done, {save_count:,} pathway saved.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pf', dest='pathway_file',
                        default='../data_1080/pathway_rc/Gene_RestrictedConcept.pathway.107.tsv',
                        help='default: ../data_1080/pathway_rc/Gene_RestrictedConcept.pathway.107.tsv')
    parser.add_argument('-hp', dest='hypothesis_path',
                        default='../data_1080/pathway_rc/hypothesis_module_for_pathway107',
                        help='default: ../data_1080/pathway_rc/hypothesis_module_for_pathway107')
    parser.add_argument('-mp', dest='match_prefix',
                        choices={'MAPT', 'APP'},
                        default='APP')
    parser.add_argument('-sf', dest='save_file',
                        required=True)
    args = parser.parse_args()

    main(args.pathway_file, args.hypothesis_path, args.match_prefix, args.save_file)
