# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 11/11/2022 16:34
@Author: yao
"""

"""
template code for db file re-name.
*.iASiS.DB.tsv
*.Pathway.DB.tsv
"""

import os
import gzip


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


base_path_dir = '/home/xzyao/AD-PNRLE/result/complete_database'

db_path_list = ['drug_db_dir', 'stage_db_dir', 'subtype_db_dir', 'db_dir']

gene_info = '/home/xzyao/AD-PNRLE/data/ncbi-data/Homo_sapiens.gene_info.gz'

id_to_symbol = read_gene_info(gene_info)

for db_dir in db_path_list:

    print(f'Process: {db_dir}')
    db_path = f'{base_path_dir}/{db_dir}'

    pathway_file = f'{db_path}/iASiS-PNRLE.Pathway.DB.tsv'
    pathway_save_file = f'{db_path}/iASiS-PNRLE.Pathway.DB.norm.tsv'

    complete_gene_file = f'{db_path}/iASiS-PNRLE.CompleteGene.DB.tsv'
    complete_save_file = f'{db_path}/iASiS-PNRLE.CompleteGene.DB.norm.tsv'


    # Complete file
    if os.path.exists(complete_gene_file):
        with open(complete_gene_file) as f, open(complete_save_file, 'w') as wf:
            head = f.readline()
            wf.write(head)

            for line in f:
                l = line.strip().split('\t')
                gene_token = l[2]
                symbol, _, entrez = gene_token.split('|')

                if id_to_symbol.get(entrez):
                    if id_to_symbol[entrez] != symbol:
                        # print(l)
                        gene_token = f'{id_to_symbol[entrez]}|Gene|{entrez}'
                        l[2] = gene_token
                        # print(l)
                        # input()
                l = '\t'.join(l)
                wf.write(f'{l}\n')


    # Pathway file
    with open(pathway_file) as f, open(pathway_save_file, 'w') as wf:

        head = f.readline()
        wf.write(head)

        for line in f:
            l = line.strip().split('\t')

            entrez = l[2]
            symbol = l[3]
            gene_token = l[6]

            if id_to_symbol.get(entrez):
                if id_to_symbol[entrez] != symbol:
                    # print(l)
                    symbol = id_to_symbol[entrez]
                    l[3] = symbol
                    l[6] = f'{symbol}|Gene|{entrez}'
                # print(l)
                # input()
            l = '\t'.join(l)

            wf.write(f'{l}\n')

            """6       3.5     2896    GRN     GO:0051402      Apoptosis of neurons    GRN|Gene|2896   ThemeOf"""




