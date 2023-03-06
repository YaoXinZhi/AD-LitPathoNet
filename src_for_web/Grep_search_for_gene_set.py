# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 17/06/2022 9:29
@Author: XINZHI YAO
"""

import os
import argparse
from tqdm import tqdm
from collections import defaultdict

def read_gene_set(gene_set_file: str):

    gene_set = set()
    with open(gene_set_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            entrez = l[0]
            gene_set.add(entrez)
    print(f'Gene set: {len(gene_set)}')
    return gene_set


def gene_set_search(gene_set_file: str, chain_file: str, save_file: str):

    gene_set = read_gene_set(gene_set_file)

    for entrez in tqdm(gene_set):

        go_commend = f'cat {chain_file}| grep -E "Gene-{entrez}.*ThemeOf.*Var.*rs[0-9]+.*CauseOf.*GO:" >> {save_file}'
        hpo_commend = f'cat {chain_file}| grep -E "Gene-{entrez}.*ThemeOf.*Var.*rs[0-9]+.*CauseOf.*HP:" >> {save_file}'
        dis_commend = f'cat {chain_file}| grep -E "Gene-{entrez}.*ThemeOf.*Var.*rs[0-9]+.*CauseOf.*MESH:" >> {save_file}'
        # print(go_commend)
        os.system(go_commend)
        os.system(hpo_commend)
        os.system(dis_commend)

    print(f'{save_file} save done.')


def main():
    parser = argparse.ArgumentParser(description='Grep search for gene set.')
    parser.add_argument('-gf', dest='gene_set_file', required=True)
    parser.add_argument('-cf', dest='chain_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    gene_set_search(args.gene_set_file, args.chain_file, args.save_file)



if __name__ == '__main__':
    main()



