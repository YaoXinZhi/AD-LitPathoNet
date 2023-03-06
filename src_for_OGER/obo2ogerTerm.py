# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 16/06/2022 9:33
@Author: XINZHI YAO
"""

import os
import re

import argparse

def obo2oger(obo_file: str, save_file: str):

    term_count = 1
    _id = ''
    name = ''
    namespace = ''
    syn_set = set()
    with open(obo_file) as f, open(save_file, 'w') as wf:
        for line in f:
            l = line.strip()
            if l.startswith('[Term]'):
                if not _id:
                    continue
                if syn_set:
                    syn_wf = ', '.join(syn_set)
                    wf.write(f'c{term_count}\tGene Ontology\t'
                             f'{_id}\t{name}\t{syn_wf}\t{namespace}\n')
                    syn_set = set()
                else:
                    wf.write(f'c{term_count}\tGene Ontology\t'
                             f'{_id}\t{name}\t{name}\t{namespace}\n')
                term_count += 1
            if l.startswith('id:'):
                _id = l.split()[1]
            if l.startswith('name:'):
                name = ' '.join(l.split()[1:])
            if l.startswith('namespace:'):
                namespace = l.split()[1]
            if l.startswith('synonym:'):
                if not 'EXACT' in l:
                    continue
                syn_set.update(re.findall(r'"(.*)"', l))
    print(f'{save_file} save done.')


def main():

    parser = argparse.ArgumentParser('obo 2 oger term file.')
    parser.add_argument('-of', dest='obo_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    obo2oger(args.obo_file, args.save_file)


if __name__ == '__main__':
    main()




