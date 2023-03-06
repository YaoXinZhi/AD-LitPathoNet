# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 02/06/2022 9:31
@Author: XINZHI YAO
"""

import os
import gzip

from collections import defaultdict

def read_gene_file(gene_file: str):

    entrez_set = set()
    with open(gene_file) as f:
        for line in f:
            l = line.strip().split('\t')
            if len(l) == 2:
                entrez = l[1]
            else:
                entrez = l[0]

            entrez_set.add(entrez)
    return entrez_set

def read_gene_to_go(gene2go_file: str):

    gene_to_go = defaultdict(set)
    with gzip.open(gene2go_file) as f:
        f.readline()
        for line in f:
            l = line.decode('utf-8').strip().split('\t')

            tax_id = l[0]
            if tax_id != '9606':
                continue

            gene_id = l[1]
            go_id = l[2]

            gene_to_go[gene_id].add(go_id)
    return gene_to_go


def read_obo(obo_file: str):
    id_to_name = {}

    # for MF, BP, CC
    id_to_type = {}
    _id = ''
    alt_id_set = set()
    no_annotation_id_set = set()
    with open(obo_file) as f:
        for line in f:

            l = line.strip().split(': ')
            if len(l) < 2:
                continue

            if line.startswith('id:'):
                _id = l[ 1 ]
                alt_id_set = set()

            if line.startswith('name:'):
                name = l[ 1 ]
                id_to_name[ _id ] = name

            if line.startswith('comment:'):
                if 'not be used for direct annotation' in line:
                    no_annotation_id_set.add(_id)
                    no_annotation_id_set.update(alt_id_set)

            if line.startswith('subset:'):
                if 'gocheck_do_not_annotate' in line:
                    no_annotation_id_set.add(_id)
                    no_annotation_id_set.update(alt_id_set)

            if line.startswith('namespace:'):
                term_type = l[ 1 ]
                id_to_type[ _id ] = term_type

            if line.startswith('alt_id:'):
                alt_id = l[ 1 ]
                alt_id_set.add(alt_id)
                id_to_name[ alt_id ] = name

    print(f'Term: {len(id_to_name):,}')
    return id_to_name, no_annotation_id_set, id_to_type

def read_gene_info(gene_info_file: str):
    entrez_to_symbol = {}
    with gzip.open(gene_info_file) as f:
        f.readline()
        for line in f:
            l = line.decode('utf-8').strip().split('\t')
            entrez = l[1]
            symbol = l[2]

            entrez_to_symbol[entrez] = symbol

    print(f'gene count: {len(entrez_to_symbol):,}')
    return entrez_to_symbol

def generate_table_html(gene_file: str):

    selected_gene_to_go, entrez_to_symbol, go_id_to_name, selected_gene_to_hpo, hpo_id_to_name = get_gene_related_go(gene_file)

    tr_element = ''
    for entrez, go_set in selected_gene_to_go.items():
        for go_id in go_set:
            go_name = go_id_to_name[go_id]
            symbol = entrez_to_symbol[entrez]
            tr_element += f"""
                            <tr>
                                <td>{entrez}</td>
                                <td>{symbol}</td>
                                <td>{go_id}</td>
                                <td>{go_name}</td>
                            </tr>
                        """

    for entrez, hpo_set in selected_gene_to_hpo.items():
        for hpo_id in hpo_set:
            hpo_name = hpo_id_to_name[hpo_id]
            symbol = entrez_to_symbol[entrez]
            tr_element += f"""
                            <tr>
                                <td>{entrez}</td>
                                <td>{symbol}</td>
                                <td>{hpo_id}</td>
                                <td>{hpo_name}</td>
                            </tr>
                        """
    return tr_element


def read_gene_to_hpo(gene_to_hpo_file: str, selected_entrez_set: set):

    selected_gene_to_hpo = defaultdict(set)
    hpo_id_to_name = {}
    with open(gene_to_hpo_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')

            hpo_id = l[0]
            hpo_name = l[1]
            entrez = l[2]
            symbol = l[3]

            if entrez not in selected_entrez_set:
                continue

            selected_gene_to_hpo[entrez].add(hpo_id)
            hpo_id_to_name[hpo_id] = hpo_name

    return selected_gene_to_hpo, hpo_id_to_name




def get_gene_related_go(gene_file: str,
                        go_file: str='../data/obo-data/go.obo',
                        gene_to_go_file: str='../data/ncbi-data/gene2go.gz',
                        gene_to_hpo_file: str='../data/ncbi-data/phenotype_to_genes.txt',
                        homo_gene_info_file: str='../data/ncbi-data/Homo_sapiens.gene_info.gz'):

    entrez_set = read_gene_file(gene_file)
    go_id_to_name, _, go_to_type = read_obo(go_file)
    entrez_to_symbol = read_gene_info(homo_gene_info_file)

    gene_to_go = read_gene_to_go(gene_to_go_file)

    selected_gene_to_hpo, hpo_id_to_name = read_gene_to_hpo(gene_to_hpo_file, entrez_set)

    selected_gene_to_go = {entrez: go_set for entrez, go_set in gene_to_go.items() if entrez in entrez_set}


    return selected_gene_to_go, entrez_to_symbol, go_id_to_name, selected_gene_to_hpo, hpo_id_to_name


if __name__ == '__main__':

    pass

