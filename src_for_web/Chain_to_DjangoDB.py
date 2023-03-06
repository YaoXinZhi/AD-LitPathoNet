# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 27/06/2022 20:29
@Author: XINZHI YAO
"""

"""
改代码用于代替原有的bash命令行
对链条结果文件进行处理
同时存成Django数据库的格式方便导入

当前版本：
Entrez Symbol MutationType MutationSubType Chain_line

候选：
Entrez Symbol MutationType MutationSubType GOID GOName HPOID HPOName DisID DisName Chain_line
1. 链条筛选 
    除了Reg都标准化
"""

import argparse
import os
import re
from collections import defaultdict
from itertools import product

class evi_class:
    def __init__(self, id_to_element: dict, chain_line: str):
        self.id_to_element = id_to_element
        self.chain_line = chain_line

def chain_parser(chain: list):
    """
    当前版本不单独存 GO HPO Dis
    1. 除了reg都必须标准化
    2. 判断是否有ThemeOf的基因
    """

    mutation_type = ''
    sub_type = ''
    gene_set = set()
    hpo_set = set()
    go_set = set()
    dis_set = set()
    reg_type = ''

    saw_reg = False

    source_gene_set = set()
    target_gene_set = set()

    rs_id = 'None'

    # entrez/bio-conceptId to element
    id_to_element = {}
    # get the key of Database
    # Gene, MutationType, MutationSubType,
    # 7-4: ThemeGene
    reg_set = {'Reg', 'NegReg', 'PosReg'}
    for idx, element in enumerate(chain):

        if element in {'ThemeOf', 'CauseOf'}:
            continue

        element = eval(element)

        if element[ 1 ] in reg_set:
            reg_type = element[ 1 ]

            # print(element)
            id_to_element[reg_type] = element
            saw_reg = True

        if element[ 1 ] not in reg_set and element[ -1 ] == '-':
            return '', '', '', '', '', '', '', '', '', '', '', rs_id, ''

        # source gene
        if element[ 1 ] == 'Gene' and not saw_reg:
            norm_set = element[ -1 ].split(', ')
            for norm in norm_set:
                symbol, entrez = norm.split('-Gene-')
                gene_set.add((entrez, symbol))
                source_gene_set.add((entrez, symbol))
                id_to_element[ entrez ] = element


        # mutation
        elif element[ 1 ] == 'Var':
            if len(element[ -1 ].split(', ')) > 1:
                return '', '', '', '', '', '', '', '', '', '', '', rs_id, ''
            mutation_type = '-'.join(element[ -1 ].split('-')[ :-2 ])
            sub_type = element[ -1 ].split('-')[ -1 ]
            rs_id = re.findall(r'rs[0-9]+', sub_type)
            if re.findall(r'rs[0-9]+', sub_type):
                mutation_type = 'point mutations'
                sub_type = 'RsMutation'
                rs_id = rs_id[ 0 ]
                id_to_element[ rs_id ] = element

            else:
                rs_id = 'None'
                id_to_element[ sub_type ] = element

        # biological concept
        else:
            norm_set = element[ -1 ].split(', ')

            for norm in norm_set:
                norm_id = norm.split('-')[ -1 ]
                norm_name = '-'.join(norm.split('-')[ :-2 ])
                id_to_element[ norm_id ] = element

                if norm_id.startswith('GO:'):
                    go_set.add((norm_id, norm_name))
                elif norm_id.startswith('HP:'):
                    hpo_set.add((norm_id, norm_name))
                elif norm_id.startswith('MESH:'):
                    dis_set.add((norm_id, norm_name))

    # 1. 判断该链条是否有ThemeOf 的基因
    have_theme_gene = 'False'
    if chain[ -2 ] == 'ThemeOf':
        last_element = eval(chain[ -1 ])
        # biological process
        bp_element = eval(chain[ -3 ])

        last_type = last_element[ 1 ]

        bp_type = bp_element[ 1 ]

        if last_type == 'Gene' and bp_type not in reg_set and last_element[ -1 ] != '-':

            have_theme_gene = 'True'
            for norm in last_element[ -1 ].split(', '):
                symbol, entrez = norm.split('-Gene-')

                target_gene_set.add((entrez, symbol))
                id_to_element[entrez] = last_element

    # source_gene == target_gene
    if source_gene_set & target_gene_set:
        self_interaction = 'True'
    else:
        self_interaction = 'False'

    return gene_set, mutation_type, sub_type, have_theme_gene, go_set, hpo_set, dis_set, \
           source_gene_set, target_gene_set, self_interaction, reg_type, rs_id, id_to_element


# fixme: only reg edge have all rich evi
def sent_to_rich(sent: str, chain: set, unique_edge: tuple, id_to_element: dict):
    type_to_color = {
        'Gene': '#ffa235',
        'Protein': '#ffa235',
        'Var': '#8932a5',
        'GO': '#90ed7d',
        'HPO': '#70ad47',
        'Disease': '#82B0D2',

        'Reg': '#DCDCDC',
        'PosReg': '#FF6666',
        'NegReg': '#00a8e1'
    }

    # print(unique_edge)
    # # print(chain)
    # print(id_to_element)
    # input()

    offset_to_color = {}
    for _id in unique_edge:
        if _id == 'None':
            continue

        element = id_to_element[_id]

        term, _type, offset, _ = element

        if offset == 'dbSNP':
            continue
        offset = eval(offset)

        start, end = offset
        if sent[start: end] != term:
            print('wrong offset')
            print(sent)
            print(offset)
            print(sent[start: end])
            print(term)
            input()

        if _id.startswith('GO:'):
            term_type = 'GO'
        elif _id.startswith('HP:'):
            term_type = 'HPO'
        elif _id.startswith('MESH'):
            term_type = 'Disease'
        else:
            term_type = _type

        offset_to_color[offset] = type_to_color[term_type]

    # print(offset_to_color)
    # input()

    sorted_offset = sorted(offset_to_color.keys(), key=lambda x: x[0], reverse=True)
    character_list = [s for s in sent]

    # start from last offset
    for offset in sorted_offset:
        start, end = offset
        color = offset_to_color[offset]

        character_list.insert(end, '</b></font>')
        character_list.insert(start, f'<font color="white" style="background:{color}"><b>')

    rich_text = ''.join(character_list)

    rich_text = rich_text.replace('-->', '--\>')

    return rich_text


def read_no_norm_mutation_file(no_norm_mutation_file: str):
    mutation_type_set = set()
    sub_mutation_type_set = set()
    with open(no_norm_mutation_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            if len(l) <5:
                continue

            mutation_type = l[1]

            sub_mutation_type = l[3]

            syno_name_set = set(l[4].split('|'))

            mutation_type_set.add(mutation_type)
            sub_mutation_type_set.add(sub_mutation_type)
            sub_mutation_type_set.update(syno_name_set)


    return mutation_type_set, sub_mutation_type_set


def unique_edge_to_search_key(unique_edge: tuple):
    # ('54664', 'Mutations', 'Reg', 'HP:0002180', 'None')

    source_entrez, mutation, relation, term_id, target_entrez = unique_edge

    search_key = f'{source_entrez}-{mutation}-{relation}-{term_id}-{target_entrez}'

    return search_key

def main(chain_file: str, save_path: str, prefix: str, no_norm_mutation_file: str):
    if not os.path.exists(save_path):
        os.mkdir(save_path)

    # mutation_type_set, sub_mutation_set = read_no_norm_mutation_file(no_norm_mutation_file)
    #
    # print(mutation_type_set)
    # print(sub_mutation_set)

    main_save_file = f'{save_path}/{prefix}.DjangoDB.tsv'

    go_save_file = f'{save_path}/{prefix}.GO.tsv'
    hpo_save_file = f'{save_path}/{prefix}.HPO.tsv'
    dis_save_file = f'{save_path}/{prefix}.DIS.tsv'

    ppi_save_file = f'{save_path}/{prefix}.PPI.tsv'

    # 20220715 rich evidence table
    evidence_save_file = f'{save_path}/{prefix}.evi.tsv'

    # unique_edge: entrez, mutation, relation, BiologicalConceptId, ThemeEntrez
    unique_edge_to_evi = defaultdict(set)
    unique_edge_to_eviid = {}
    # unique_edge_to_id_to_elemet = {}

    saved_ppi_set = set()
    with open(chain_file) as f, open(main_save_file, 'w') as wf, \
            open(go_save_file, 'w') as wf_go, \
            open(hpo_save_file, 'w') as wf_hpo, \
            open(dis_save_file, 'w') as wf_dis, \
            open(ppi_save_file, 'w') as wf_ppi:
        wf.write(f'ProteinOneEntrez\tProteinOneSymbol\t'
                 f'ProteinTwoEntrez\tProteinTwoSymbol\t'
                 f'MutationType\tMutationSubType\t'
                 f'HaveThemeGene\t'
                 f'SelfInteraction\t'
                 f'RegType\t'
                 f'Chain_line\n')

        wf_go.write(f'GO ID\tGO Term\t'
                    f'MutationType\tMutationSubType\t'
                    f'HaveThemeGene\t'
                    f'SelfInteraction\t'
                    f'RegType\t'
                    f'Chain_line\n')

        wf_hpo.write(f'HPO ID\tHPO Term\t'
                     f'MutationType\tMutationSubType\t'
                     f'HaveThemeGene\t'
                     f'SelfInteraction\t'
                     f'RegType\t'
                     f'Chain_line\n')

        wf_dis.write(f'Mesh ID\tMesh Term\t'
                     f'MutationType\tMutationSubType\t'
                     f'HaveThemeGene\t'
                     f'SelfInteraction\t'
                     f'RegType\t'
                     f'Chain_line\n')

        wf_ppi.write(f'Protein1 Entrez\tProtein1 Symbol\t'
                     f'Protein2 Entrez\tProtein2 Symbol\t'
                     f'SelfInteraction\t'
                     f'RegType\t'
                     f'MutationType\tMutationSubType\t'
                     f'Chain_line\n')

        for line in f:
            l = line.strip().split('\t')
            if l == [ '' ]:
                continue

            # pmid = l[ 0 ]
            # sent_id = l[ 1 ]
            # sentence = l[ 2 ]
            chain = l[ 3: ]

            gene_set, mutation_type, sub_type, \
            have_theme_gene, \
            go_set, hpo_set, dis_set, \
            source_gene_set, target_gene_set, \
            self_interaction, \
            reg_type, rs_id, \
            id_to_element = chain_parser(chain)

            # if (mutation_type not in mutation_type_set or sub_type not in sub_mutation_set) and sub_type!='RsMutation':
            #     # print(f'mutation: {mutation_type}')
            #     # print(f'sub mutation: {sub_type}')
            #     # input()
            #     mutation_type = 'point mutation'
            #     sub_type = 'Mutations'


            if not gene_set \
                    or not mutation_type \
                    or not sub_type:
                continue

            if have_theme_gene == 'False':
                for (entrez, symbol) in gene_set:
                    wf.write(f'{entrez}\t{symbol}\t'
                             f'None\tNone\t'
                             f'{mutation_type}\t{sub_type}\t'
                             f'{have_theme_gene}\t'
                             f'{self_interaction}\t'
                             f'{reg_type}\t'
                             f'{line.strip()}\n')

            for (go_id, go_term) in go_set:
                wf_go.write(f'{go_id}\t{go_term}\t'
                            f'{mutation_type}\t{sub_type}\t'
                            f'{have_theme_gene}\t'
                            f'{self_interaction}\t'
                            f'{reg_type}\t'
                            f'{line.strip()}\n')

            for (hpo_id, hpo_term) in hpo_set:
                wf_hpo.write(f'{hpo_id}\t{hpo_term}\t'
                             f'{mutation_type}\t{sub_type}\t'
                             f'{have_theme_gene}\t'
                             f'{self_interaction}\t'
                             f'{reg_type}\t'
                             f'{line.strip()}\n')

            for (dis_id, dis_term) in dis_set:
                wf_dis.write(f'{dis_id}\t{dis_term}\t'
                             f'{mutation_type}\t{sub_type}\t'
                             f'{have_theme_gene}\t'
                             f'{self_interaction}\t'
                             f'{reg_type}\t'
                             f'{line.strip()}\n')

            if have_theme_gene == 'True' and source_gene_set and target_gene_set:
                for source_gene, target_gene in product(source_gene_set, target_gene_set):

                    if (source_gene, target_gene) in saved_ppi_set \
                            or (target_gene, source_gene) in saved_ppi_set:
                        continue
                    saved_ppi_set.add((source_gene, target_gene))

                    source_entrez, source_symbol = source_gene
                    target_entrez, target_symbol = target_gene

                    wf.write(f'{source_entrez}\t{source_symbol}\t'
                             f'{target_entrez}\t{target_symbol}\t'
                             f'{mutation_type}\t{sub_type}\t'
                             f'{have_theme_gene}\t'
                             f'{self_interaction}\t'
                             f'{reg_type}\t'
                             f'{line.strip()}\n')

                    wf_ppi.write(f'{source_entrez}\t{source_symbol}\t'
                                 f'{target_entrez}\t{target_symbol}\t'
                                 f'{self_interaction}\t'
                                 f'{reg_type}\t'
                                 f'{mutation_type}\t{sub_type}\t'
                                 f'{line.strip()}\n')

            # add unique edge data
            if not target_gene_set:
                target_gene_set = {('None', 'None')}

            for (source_entrez, source_symbol) in source_gene_set:
                for (target_entrez, target_symbol) in target_gene_set:

                    for (go_id, go_term) in go_set:
                        if rs_id != 'None':
                            unique_edge = (source_entrez, rs_id, reg_type, go_id, target_entrez)
                        else:
                            unique_edge = (source_entrez, sub_type, reg_type, go_id, target_entrez)

                        unique_edge_to_evi[ unique_edge ].add(evi_class(id_to_element, line.strip()))
                        if not unique_edge_to_eviid.get(unique_edge):
                            unique_edge_to_eviid[ unique_edge ] = len(unique_edge_to_eviid)

                    for (hpo_id, hpo_term) in hpo_set:
                        if rs_id != 'None':
                            unique_edge = (source_entrez, rs_id, reg_type, hpo_id, target_entrez)
                        else:
                            unique_edge = (source_entrez, sub_type, reg_type, hpo_id, target_entrez)

                        unique_edge_to_evi[ unique_edge ].add(evi_class(id_to_element, line.strip()))

                        if not unique_edge_to_eviid.get(unique_edge):
                            unique_edge_to_eviid[ unique_edge ] = len(unique_edge_to_eviid)

                    for (dis_id, dis_term) in dis_set:
                        if rs_id != 'None':
                            unique_edge = (source_entrez, rs_id, reg_type, dis_id, target_entrez)
                        else:
                            unique_edge = (source_entrez, sub_type, reg_type, dis_id, target_entrez)

                        unique_edge_to_evi[ unique_edge ].add(evi_class(id_to_element, line.strip()))

                        if not unique_edge_to_eviid.get(unique_edge):
                            unique_edge_to_eviid[ unique_edge ] = len(unique_edge_to_eviid)


    with open(evidence_save_file, 'w') as wf_evi:
        # evidence_id not unique
        wf_evi.write('SearchID\t'
                     'UniqueEdge\t'
                     'PMID\tSentID\t'
                     'Sentence\t'
                     'RichEvidence\n')
        for unique_edge, evi_set in unique_edge_to_evi.items():

            saved_evi = set()
            search_key = unique_edge_to_search_key(unique_edge)
            for evi_element in evi_set:

                evi = evi_element.chain_line
                id_to_element = evi_element.id_to_element

                l = evi.split('\t')

                pmid = l[ 0 ]
                sent_id = l[ 1 ]
                sentence = l[ 2 ]
                chain = l[ 3: ]

                if (pmid, sentence) in saved_evi:
                    continue
                saved_evi.add((pmid, sentence))

                rich_sent = sent_to_rich(sentence, chain, unique_edge, id_to_element)

                wf_evi.write(f'{search_key}\t{unique_edge}\t{pmid}\t{sent_id}\t{sentence}\t{rich_sent}\n')

    print(f'{main_save_file} save done.')
    print(f'{go_save_file} save done.')
    print(f'{hpo_save_file} save done.')
    print(f'{dis_save_file} save done.')
    print(f'{ppi_save_file} save done.')
    print(f'{evidence_save_file} save done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Chain to Django DB format.')
    parser.add_argument('-if', dest='input_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    parser.add_argument('-pr', dest='prefix', required=True)

    parser.add_argument('-nv',dest='no_norm_var_file', default='/home/xzyao/AD-PNRLE/data/Non-Normalized_Mutations/MutationTypeDictionary20220711.txt',
                        help='/home/xzyao/AD-PNRLE/data/Non-Normalized_Mutations/MutationTypeDictionary20220711.txt')
    args = parser.parse_args()

    main(args.input_file, args.save_path, args.prefix, args.no_norm_var_file)
