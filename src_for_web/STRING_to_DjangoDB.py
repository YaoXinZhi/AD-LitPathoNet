# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 04/07/2022 20:53
@Author: XINZHI YAO
"""

import os
import gzip
import argparse
from collections import defaultdict

def read_gene_info(gene_info_file: str):

    entrez_to_symbol = {}
    with gzip.open(gene_info_file) as f:
        f.readline()
        for line in f:
            l = line.decode('utf-8').strip().split('\t')

            if l[0] != '9606':
                continue
            entrez = l[1]
            symbol = l[2]

            entrez_to_symbol[entrez] = symbol
    return entrez_to_symbol

def read_ensembl_to_entrez(ensembl_to_entrez_file: str):

    ensembl_to_entrez = {}
    with gzip.open(ensembl_to_entrez_file) as f:
        f.readline()
        for line in f:
            l = line.decode('utf-8').strip().split('\t')

            if l[0] != '9606':
                continue

            ensembl = l[-1].split('.')[0]
            entrez = l[1]

            ensembl_to_entrez[ensembl] = entrez
    return ensembl_to_entrez


def main(db_file: str, save_file: str,
         ensembl_to_entrez_file: str, entrez_to_symbol_file: str):

    ensembl_to_entrez = read_ensembl_to_entrez(ensembl_to_entrez_file)

    entrez_to_symbol = read_gene_info(entrez_to_symbol_file)

    save_count = 0
    miss_count = 0
    with gzip.open(db_file) as f, open(save_file, 'w') as wf:
        wf.write(f'protein1_entrez\tprotein1_symbol\t'
                 f'protein2_entrez\tprotein2_symbol\t'
                 f'combined_score\n')
        l = f.readline().decode('utf-8').strip().split()
        # print(l)
        # print(l[9], l[11])
        for line in f:
            l = line.decode('utf-8').strip().split()
            # 9606.ENSP00000000233
            protein_1_ensembl = l[0].split('.')[1]
            protein_2_ensembl = l[1].split('.')[1]

            experiment = int(l[9])
            database = int(l[11])

            score = l[-1]

            if experiment == 0 or database == 0:
                continue

            if ensembl_to_entrez.get(protein_1_ensembl) and ensembl_to_entrez.get(protein_2_ensembl):
                entrez_1 = ensembl_to_entrez[protein_1_ensembl]
                entrez_2 = ensembl_to_entrez[protein_2_ensembl]

                symbol_1 = entrez_to_symbol[entrez_1]
                symbol_2 = entrez_to_symbol[entrez_2]
            else:
                miss_count += 1
                continue

            wf.write(f'{entrez_1}\t{symbol_1}\t'
                     f'{entrez_2}\t{symbol_2}\t'
                     f'{score}\n')
            save_count += 1

    print(f'{save_count} PPI saved, {miss_count} missed')
    print(f'{save_file} save done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='STRING Filter for Django DB')
    parser.add_argument('-df', dest='string_db_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    parser.add_argument('-ef', dest='ensembl_to_entrez_file',
                        default='/mnt/disk2/xzyao/ncbi_data/gene2ensembl.gz')
    parser.add_argument('-yf', dest='entrez_to_symbol_file',
                        default='/mnt/disk2/xzyao/ncbi_data/gene_info.gz')
    args = parser.parse_args()

    main(args.string_db_file, args.save_file,
         args.ensembl_to_entrez_file, args.entrez_to_symbol_file)


