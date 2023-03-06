# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 20/09/2022 16:09
@Author: XINZHI YAO
"""

import os
from tqdm import tqdm

from collections import defaultdict

def go_umls_mapping(umls_file: str, go_mapping_save_file: str,
                    hpo_mapping_save_file: str,
                    mesh_umls_mapping_save_file: str):
    """
    GO条目和UMLS条目存在多对多
    """

    umls_set = set()
    total_go_set = set()
    total_hpo_set = set()
    total_mesh_set = set()

    umls_to_go = defaultdict(set)
    umls_to_hpo = defaultdict(set)
    umls_to_mesh = defaultdict(set)
    with open(umls_file) as f:
        for line in tqdm(f):
            l = line.strip().split('|')
            cui = l[0]
            alt_id = l[10]
            if alt_id.startswith('GO:'):

                umls_to_go[cui].add(alt_id)
                umls_set.add(cui)
                total_go_set.add(alt_id)

            if alt_id.startswith('HP:'):
                umls_to_hpo[cui].add(alt_id)
                umls_set.add(cui)
                total_hpo_set.add(alt_id)

            if l[11] == 'MSH' and alt_id.startswith('D'):
                umls_to_mesh[cui].add(alt_id)
                umls_set.add(cui)
                total_mesh_set.add(alt_id)

    with open(go_mapping_save_file, 'w') as wf:
        wf.write(f'CUI\tGO\n')
        for umls, go_set in umls_to_go.items():
            for go in go_set:
                wf.write(f'{umls}\t{go}\n')

    with open(hpo_mapping_save_file, 'w') as wf:
        wf.write(f'CUI\tHPO\n')
        for umls, hpo_set in umls_to_hpo.items():
            for hpo in hpo_set:
                wf.write(f'{umls}\t{hpo}\n')

    with open(mesh_umls_mapping_save_file, 'w') as wf:
        wf.write(f'CUI\tMeSH\n')
        for umls, mesh_set in umls_to_mesh.items():
            for mesh in mesh_set:
                wf.write(f'{umls}\t{mesh}\n')

    print(f'"{go_mapping_save_file}" and "{hpo_mapping_save_file}" saved, '
          f'includes {len(total_go_set):,} GO terms, '
          f'{len(total_hpo_set)} HPO terms,'
          f'{len(total_mesh_set)} MeSH terms and {len(umls_set):,} UMLS terms.')


def main():

    umls_file = '/home/xzyao/AD-PNRLE/data/umls/MRCONSO.RRF'

    # 1. go-umls mapping
    go_file = '/home/xzyao/AD-PNRLE/data/obo-data/go.obo'
    go_umls_mapping_save_file = '/home/xzyao/AD-PNRLE/data/go_hpo_umls_mapping/go_umls_mapping.tsv'


    # 2. hpo-umls mapping
    hpo_file = '/home/xzyao/AD-PNRLE/data/obo-data/hp.obo'
    hpo_umls_mapping_save_file = '/home/xzyao/AD-PNRLE/data/go_hpo_umls_mapping/hpo_umls_mapping.tsv'

    # 3. mesh_mapping
    mesh_umls_mapping_save_file = '/home/xzyao/AD-PNRLE/data/go_hpo_umls_mapping/mesh_umls_mapping.tsv'

    go_umls_mapping(umls_file, go_umls_mapping_save_file, hpo_umls_mapping_save_file, mesh_umls_mapping_save_file)



if __name__ == '__main__':
    main()