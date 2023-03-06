# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 20/11/2022 22:27
@Author: yao
"""
import os.path

"""
改代码将合并后的GO HPO Dis数据库合并成一个
并且添加 Entrez Symbol列
"""


def main(db_path: str):

    go_file = f'{db_path}/iASiS-PNRLE.GO.tsv'
    hpo_file = f'{db_path}/iASiS-PNRLE.HPO.tsv'
    dis_file = f'{db_path}/iASiS-PNRLE.DIS.tsv'

    save_file = f'{db_path}/PNRLE.BioConcept.tsv'

    saved_set = set()
    with open(save_file, 'w') as wf:
        wf.write(f'Entre\tSymbol\t'
                 f'TermID\tTerm\t'
                 f'MutationType\tMutationSubType\t'
                 f'RegType\tChain_line\n')
        for _file in [go_file, hpo_file, dis_file]:
            print(f'Processing: {os.path.basename(_file)}.')
            with open(_file) as f:
                f.readline()
                for line in f:
                    l = line.strip().split('\t')

                    term_id = l[0]
                    term = l[1]
                    mutation_type = l[2]
                    mutation_sub_type = l[3]
                    reg = l[6]
                    pmid = l[7]
                    sent = l[9]

                    chain = '\t'.join(l[7:])

                    gene = eval(l[10])[-1]
                    #print(gene)
                    entrez = gene.split('-')[-1]
                    symbol = '-'.join(gene.split('-')[:-2])

                    uniq_chain = (entrez, term_id, pmid, sent)
                    if uniq_chain in saved_set:
                        continue
                    saved_set.add(uniq_chain)

                    wf.write(f'{entrez}\t{symbol}\t'
                             f'{term_id}\t{term}\t'
                             f'{mutation_type}\t{mutation_sub_type}\t'
                             f'{reg}\t{chain}\n')

    print(f'{os.path.basename(save_file)} save done, include {len(saved_set):,} chains.')


    # iASiS - PNRLE.DIS.tsv
    # iASiS - PNRLE.DjangoDB.tsv
    # iASiS - PNRLE.GO.tsv
    # iASiS - PNRLE.HPO.tsv



if __name__ == '__main__':
    merge_db_path = '/Users/yao/Nutstore Files/Mac2PC/AD-PNRLE/result/complete_database/base_db_merge_dir'

    main(merge_db_path)



