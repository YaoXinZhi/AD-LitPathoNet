# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 30/09/2022 16:03
@Author: XINZHI YAO
"""

import os
import argparse
import re
from itertools import product
from collections import defaultdict


def token_split(pathway_token: str):
    name, _type, _id = pathway_token.split('|')
    return name, _type, _id

def extract_triple_from_pathway(path_line: str, label_line: str, type_line: str,
                                rich_evidence: tuple, pathway_concept_set: set,
                                syno_to_mutation_type: dict, syno_to_sub_mutation_type: dict):

    # triple_set used to filter the iASiS evidence
    triple_set = set()

    path_line_split = path_line.split('\t')

    pathway_score = path_line_split[1]
    pathway_length = path_line_split[2]

    path_list = path_line_split[5:]
    label_list = label_line.split('\t')[5:]
    type_list = type_line.split('\t')[5:]

    # fine mutation_type/sub_mutation_type
    mutation_idx = type_list.index('Var')
    if mutation_idx == -1:
        print(f'path_list: {path_list}')
        print(f'label_list: {label_list}')
        print(f'type_list: {type_list}')
        raise ValueError('no mutation in the path.')
    mutation_id = label_list[mutation_idx]
    if re.findall(r'rs[0-9]+', mutation_id):
        mutation_type = 'point mutations'
        mutation_sub_type = 'RsMutation'
    else:
        mutation_syno = mutation_id.split('-')[-1]
        mutation_type = syno_to_mutation_type[mutation_syno]
        mutation_sub_type = syno_to_sub_mutation_type[mutation_syno]

    # pathway compound token
    pathway_token_list = []

    # compound token list for PNRLE part in pathway
    pnrle_token_list = []

    pmid, sent_id, rich_sent, have_dbsnp_theme = rich_evidence

    offset = 0
    while offset+2 < len(path_list):
        source = path_list[offset]
        source_id = label_list[offset]
        source_type = type_list[offset]

        relation = path_list[offset+1]

        target = path_list[offset+2]
        target_id = label_list[offset+2]
        target_type = type_list[offset+2]

        if offset == 0:
            token_compound = f'{source}|{source_type}|{source_id}'
            pathway_token_list.append(token_compound)
            pnrle_token_list.append(token_compound)

        token_compound = f'{target}|{target_type}|{target_id}'
        pathway_token_list.append(relation)
        pathway_token_list.append(token_compound)

        offset += 2
        # pnrle part edge
        if relation in {'ThemeOf', 'Reg', 'PosReg', 'NegReg'}:
            pnrle_token_list.append(relation)
            pnrle_token_list.append(token_compound)
            continue
        triple_set.add((source_id, relation, target_id))

    # start and end concept in pathway
    start_token, _, start_id = token_split(pathway_token_list[0])
    end_token, _, end_id = token_split(pathway_token_list[-1])

    # start and end concept in pnrle part
    pnrle_end_token, _, pnrle_end_id = token_split(pnrle_token_list[-1])

    if pnrle_end_id in pathway_concept_set:
        pnrle_token_wf = '\t'.join(pnrle_token_list)
        pnrle_wf = f'{pnrle_end_id}\t{pnrle_end_token}\t' \
                   f'{pnrle_token_wf}\t' \
                   f'{pmid}\t{sent_id}\t{rich_sent}\t' \
                   f'{have_dbsnp_theme}'
    else:
        pnrle_wf = ''


    pathway_token_wf = '\t'.join(pathway_token_list)
    # pathway_wf is used to write into Pathway DB.
    pathway_wf = f'{pathway_length}\t{pathway_score}\t' \
                 f'{start_id}\t{start_token}\t' \
                 f'{end_id}\t{end_token}\t' \
                 f'{pathway_token_wf}\t{pmid}\t{sent_id}\t{rich_sent}' \
                 f'\t{have_dbsnp_theme}\t' \
                 f'{mutation_type}\t{mutation_sub_type}'

    return triple_set, pathway_wf, pnrle_wf

def return_unique_edge(path_line: str, label_line: str):
    path_list = [ elem for elem in path_line.split('\t')[ 5: ] if elem != '-' ]

    label_list = [ label for label in label_line.split('\t')[ 5: ] if label != '-' ]

    relation = path_list[ 3 ]


    if relation not in {'Reg', 'PosReg', 'NegReg'}:
        print('not in relation')
        print(path_line)
        input()

    gene = label_list[0]
    var = label_list[ 1 ]
    concept = label_list[ 2 ]
    theme_gene = 'None'

    #  Non normalized variation
    var_id = var
    key_unique_edge = (gene, var_id, relation, concept, theme_gene)

    unique_edge = (gene, var, relation, concept, theme_gene)

    return unique_edge, key_unique_edge

def extract_triple_from_pathway_file(pathway_file: str, pnrle_edge_to_evidence: str,
                                     pathway_concept_set: set,
                                     syno_to_mutation_type: dict, syno_to_sub_mutation_type: dict):

    # 10-6 other related gene
    # concept_set = set()

    all_triple_set = set()
    pathway_wf_set = set()
    # pnrle part in pathway, for other related gnee
    pnrle_wf_set = set()
    miss_count = 0
    missed_concept = set()
    match_count = 0
    with open(pathway_file) as f:
        f.readline()
        while True:
            path_line = f.readline().strip()
            type_line = f.readline().strip()
            label_line = f.readline().strip()

            line = f.readline()
            if not line:
                break

            unique_edge, key_unique_edge = return_unique_edge(path_line, label_line)

            if not pnrle_edge_to_evidence.get(key_unique_edge):

                print('missed')
                # print(pnrle_edge_to_evidence.keys())
                # print(key_unique_edge)
                # print(path_line)
                # input()

                missed_concept.add(list(key_unique_edge)[-2])
                miss_count += 1
                continue

            rich_evidence = pnrle_edge_to_evidence[key_unique_edge]
            triple_set, pathway_wf, pnrle_wf = extract_triple_from_pathway(path_line, label_line, type_line,
                                                                 rich_evidence, pathway_concept_set,
                                                                           syno_to_mutation_type, syno_to_sub_mutation_type)
            all_triple_set.update(triple_set)
            pathway_wf_set.add(pathway_wf)
            if pnrle_wf:
                pnrle_wf_set.add(pnrle_wf)
            match_count += 1

    return all_triple_set, pathway_wf_set, pnrle_wf_set

def get_evi(relation: str):
    resource = re.findall(r'resource: \[\'(.*?)\']', relation)
    if not resource:
        resource = 'None'
    else:
        resource = resource[0]

    sent_id = re.findall(r'sent_id: \[\'(.*?)\']', relation)
    if not sent_id:
        sent_id = 'None'
    else:
        sent_id = sent_id[0]
    return resource, sent_id

def save_pathway_db(pathway_db_save_file: str, pathway_wf_set: set):

    with open(pathway_db_save_file, 'w') as wf:
        wf.write(f'PathLength\tScore\tStartID\tStartConcept\t'
                 f'EndID\tEndConcept\t'
                 f'Path\t'
                 f'PMID\tSentID\tRiceEvidence\t'
                 f'HavedbsnpTheme\tMutationType\tMutationSubType\n')
        for pathway in pathway_wf_set:

            # replace '".*?"' to '.*?'
            if re.findall(r'".*?"', pathway):
                for err_token in re.findall(r'".*?"', pathway):
                    pathway = pathway.replace(err_token, err_token[1: -1])

            # if '"Dementia, Vascular"' in pathway:
            #     print(f'replaced err label.')
            #     pathway = pathway.replace('"Dementia, Vascular"', "Dementia, Vascular")
            # if '"Diagnosis, Psychiatric"' in pathway:
            #     print(f'replaced err label.')
            #     pathway = pathway.replace('"Diagnosis, Psychiatric"', "Diagnosis, Psychiatric")

            wf.write(f'{pathway}\n')
    print(f'{pathway_db_save_file} save done, {len(pathway_wf_set):,} saved.')


def save_CompleteGene_db(save_file: str, pnrle_complete_gene_wf_set: set):

    with open(save_file, 'w') as wf:
        wf.write(f'CompleteID\tCompleteConcept\t'
                 f'Path\t'
                 f'PMID\tSentID\tRiceEvidence\t'
                 f'HavedbsnpTheme\n')
        for pnrle_wf in pnrle_complete_gene_wf_set:
            wf.write(f'{pnrle_wf}\n')
    print(f'{save_file} save done, {len(pnrle_complete_gene_wf_set):,} saved.')


def read_multi_id_file(multi_id_file: str):
    multi_id_mapping_dict = {}
    with open(multi_id_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            # concept = l[0]
            ids = l[1:]

            selected_id = ''
            for _id in ids:
                if _id.startswith('GO'):
                    selected_id = _id
                    continue
                elif _id.startswith('HP') and not selected_id.startswith('GO'):
                    selected_id = _id
                elif (_id.startswith('D') or _id.startswith('C')) and not (selected_id.startswith('GO') or selected_id.startswith('HP')):
                    selected_id = _id

            for _id in ids:
                multi_id_mapping_dict[_id] = selected_id
    return multi_id_mapping_dict



def chain_filter(chain: list):

    ban_type = set()
    # have theme of gene for last concept
    if len(chain) > 5:
        try:
            if chain[-2] == 'ThemeOf' and eval(chain[-1])[1]=='Gene':
                return False
        except:
            return False
        # gene-var 需要看作整体
        # 有些 go themeOf HPO/Disease的不能作为通路进一步连接的依据
        if chain[-4] != 'CauseOf':
            return False
    # only rs mutation
    # if not re.findall(r'rs\d+', str(chain)):
    #     return False

    # concept level filter
    for element in chain:
        if element in {'ThemeOf', 'CauseOf'}:
            continue
        # Disease in ban_type could HPO
        if eval(element)[1] in ban_type:
            if eval(element)[1] == 'Disease':
                if re.findall(r'MESH:.*?', element):
                    return False
                else:
                    return True
            else:
                return False

    # remove AD in pnrle
    if 'D000544' in chain[-1] \
        or 'HP:0002511' in chain[-1]:
        return False
    return True

def sent_to_rich(sent:str, element_list: list, unique_edge_list: list):
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

    offset_to_color = {}

    have_dbsnp_theme = 'false'
    for idx, element in enumerate(element_list):
        term, _type, offset, _ = element
        _id = unique_edge_list[idx]

        if offset == 'dbSNP':
            have_dbsnp_theme = 'true'
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
        elif _id.startswith('MESH') or _id.startswith('OMIM'):
            term_type = 'Disease'
        elif _id.startswith("C") or _id.startswith('D'):
            term_type = 'Disease'
        elif _type == 'Protein':
            term_type = 'Gene'
        else:
            if _type not in {'Var', 'Reg', 'PosReg', 'NegReg', 'Gene', 'Protein'}:
                print('wrong term type')
                print(sent)
                print(element_list)
                print(element)
                print(_id)
                print(_type)
                input()
            term_type = _type

        offset_to_color[offset] = type_to_color[term_type]

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

    return rich_text, have_dbsnp_theme


def chain_parser(chain_file: str, multi_id_mapping_dict: dict,
                 rs_to_entrez_symbol: dict):

    pnrle_edge_to_evidence = {}
    with open(chain_file) as f:
        for line in f:
            l = line.strip().split('\t')

            if l == ['']:
                continue

            # sent_id = l[0]
            pmid = l[0]
            sent_id = l[1]
            sent = l[2]

            chain = l[3:]

            if not chain_filter(chain):
                continue

            if len(chain) != 7:
                continue

            gene, _, var, _, reg, _, concept = chain[:7]

            # 1020 filter wrong chain
            if eval(gene)[1] != 'Gene' \
                    or eval(var)[1] != 'Var' \
                    or eval(reg)[1] not in {'Reg', 'PosReg', 'NegReg'}:
                continue


            gene_norm = eval(gene)[-1]

            relation_type = eval(reg)[1]
            concept_norm = eval(concept)[-1]

            var_id = eval(var)[-1].split('-')[-1]

            if gene_norm == '-' or var_id == '-' or concept_norm == '-'\
                    or not gene_norm or not var_id or not concept_norm:
                continue


            if var_id.lower().startswith('rs') and rs_to_entrez_symbol.get(var_id) and gene_norm=='-':
                entrez, symbol = rs_to_entrez_symbol[var_id]
                gene_norm = f'{symbol}-Gene-{entrez}'
                gene = str((symbol, 'Gene', 'dbSNP', gene_norm))


            gene_id = gene_norm.split(', ')[0].split('-')[-1]
            if multi_id_mapping_dict.get(gene_id):
                gene_id = multi_id_mapping_dict[gene_id]

            var_id = f'{gene_id}:Var-{var_id}'

            concept_id_set = set()
            for norm in concept_norm.split(', '):
                if len(norm.split('-')) <3:
                    continue
                concept_id = norm.split('-')[-1]

                if concept_id.startswith('MESH:'):
                    concept_id = concept_id[5:]
                elif concept_id.startswith('OMIM'):
                    continue

                if multi_id_mapping_dict.get(concept_id):
                    concept_id_set.add(multi_id_mapping_dict[concept_id])
                else:
                    concept_id_set.add(concept_id)

            element_list = [eval(gene), eval(var), eval(reg), eval(concept)]

            # for gene_id, concept_id in product(gene_id_set, concept_id_set):
            for concept_id in concept_id_set:

                unique_edge = [gene_id, var_id, relation_type, concept_id]

                rich_sent, have_dbsnp_theme = sent_to_rich(sent, element_list, unique_edge)

                pnrle_edge_to_evidence[(gene_id, var_id, relation_type, concept_id, 'None')] = (pmid, sent_id, rich_sent, have_dbsnp_theme)

    return pnrle_edge_to_evidence

def read_concept_count_file(concept_count_file: str):
    concept_set = set()
    with open(concept_count_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            concept_id = l[0]
            if concept_id.startswith('GO:') or concept_id.startswith('HP:')\
                or concept_id.startswith('D') or concept_id.startswith('C'):
                concept_set.add(concept_id)
    print(f'{len(concept_set):,} concepts in pathway file.')
    return concept_set

def read_rs_entrez_symbol_file(rs_entrez_symbol_file: str):
    rs_to_entrez_symbol = {}
    with open(rs_entrez_symbol_file) as f:
        f.readline()
        for line in f:
            rs, entrez, symbol = line.strip().split('\t')
            rs_to_entrez_symbol[rs] = (entrez, symbol)
    return rs_to_entrez_symbol

def read_no_norm_mutation_file(no_norm_mutation_file: str):

    syno_to_mutation_type = {}
    syno_to_sub_mutation_type = {}

    with open(no_norm_mutation_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            if len(l) <5:
                continue

            mutation_type = l[1]

            sub_mutation_type = l[3]

            syno_name_set = set(l[4].split('|'))

            for syno in syno_name_set:
                syno_to_mutation_type[syno] = mutation_type
                syno_to_sub_mutation_type[syno] = sub_mutation_type

                syno_to_mutation_type[sub_mutation_type] = mutation_type
                syno_to_sub_mutation_type[sub_mutation_type] = sub_mutation_type

    return syno_to_mutation_type, syno_to_sub_mutation_type


def main(pathway_file: str, iasis_file: str, save_path: str, convert_option: str,
         chain_file: str, multi_id_file: str, concept_count_file: str,
         rs_entrez_symbol_file: str, prefix: str, no_norm_mutation_file: str):

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    iasis_db_save_file = f'{save_path}/{prefix}.iASiS.DB.tsv'
    pathway_db_save_file = f'{save_path}/{prefix}.Pathway.DB.tsv'
    pnrle_complete_gene_file = f'{save_path}/{prefix}.CompleteGene.DB.tsv'

    pathway_concept_set = read_concept_count_file(concept_count_file)

    rs_to_entrez_symbol = read_rs_entrez_symbol_file(rs_entrez_symbol_file)

    print('Reading multi id file.')
    multi_id_mapping_dict = read_multi_id_file(multi_id_file)

    print('Reading PNRLE chain file.')
    pnrle_edge_to_evidence = chain_parser(chain_file, multi_id_mapping_dict, rs_to_entrez_symbol)

    print('Read no-norm mutation file.')
    syno_to_mutation_type, syno_to_sub_mutation_type = read_no_norm_mutation_file(no_norm_mutation_file)

    print('Extracting triple from pathway file.')
    triple_set, pathway_wf_set, pnrle_wf_set = extract_triple_from_pathway_file(pathway_file, pnrle_edge_to_evidence,
                                                                                pathway_concept_set,
                                                                                syno_to_mutation_type, syno_to_sub_mutation_type)


    print('Converting the iASiS and Pathway data to DB file.')
    if convert_option in {'iasis', 'all'}:
        print(f'Saving iasis database.')
        save_count = 0
        with open(iasis_file) as f, open(iasis_db_save_file, 'w') as wf:
            wf.write(f'Relation\tSource\tTarget\tResource\tSent_id\n')
            for line in f:
                l = line.strip().split('\t')
                source = l[-7]
                target = l[-3]
                relation_type = l[0]
                relation = l[1]

                resource, sent_id = get_evi(relation)
                if (source, relation_type, target) in triple_set:
                    save_count += 1
                    wf.write(f'{relation_type}\t{source}\t{target}\t{resource}\t{sent_id}\n')

        print(f'{iasis_db_save_file} save done, {save_count} iASiS triple saved.')

    if convert_option in {'pathway', 'all'}:
        print(f'Saving pathway database.')
        save_pathway_db(pathway_db_save_file, pathway_wf_set)

    #  save pnrle complete gene db file.
    if convert_option in {'supp', 'all'}:
        save_CompleteGene_db(pnrle_complete_gene_file, pnrle_wf_set)

if __name__ == '__main__':

    #如果直接从iASiS文件转 可能导致数据太大 所以从pathway文件转 只保留用到的
    parser = argparse.ArgumentParser()

    parser.add_argument('-pf', dest='pathway_file', required=True)

    parser.add_argument('-cf', dest='chain_file', required=True)

    parser.add_argument('-sp', dest='save_path', required=True)

    parser.add_argument('-pr', dest='prefix', default='iASiS-PNRLE',
                        help='default: iASiS-PNRLE')

    parser.add_argument('-if', dest='iasis_file',
                        default='../data/iASiS_all_convert/PNRLE_from_iASiS.all.FunRel.tsv',
                        help='default: ../data/iASiS_all_convert/PNRLE_from_iASiS.all.FunRel.tsv')

    # Generated by Concept_statistic.py
    parser.add_argument('-mi', dest='multi_id_file',
                        default='../result/multi_id_concept.tsv',
                        help='Generated by Concept_statistic.py, default: "../result/multi_id_concept.tsv"')

    # Generated by Concept_statistic.py
    parser.add_argument('-cc', dest='concept_count_file',
                        default='../result/concept_count.tsv',
                        help='For PNRLE supplementary db. Generated by Concept_statistic.py, default: "../result/concept_count.tsv"')

    parser.add_argument('-co', dest='convert_option',
                        choices=['pathway', 'iasis', 'supp', 'all'],
                        default='all')

    parser.add_argument('-rg', dest='rs_entrez_symbol_file',
                        default='../result/rs_to_entrez/rs-entrez-symbol.merge.tsv',
                        help='check the rs mutation with correct gene. default: ../result/rs_to_entrez/rs-entrez-symbol.merge.tsv')

    parser.add_argument('-nm', dest='no_norm_mutation_file',
                        default='/home/xzyao/AD-PNRLE/data/Non-Normalized_Mutations/MutationTypeDictionary20220711.txt',
                        help='/home/xzyao/AD-PNRLE/data/Non-Normalized_Mutations/MutationTypeDictionary20220711.txt')


    # 10-5 add rich evidence for each edge

    args = parser.parse_args()

    main(args.pathway_file, args.iasis_file, args.save_path, args.convert_option,
         args.chain_file, args.multi_id_file, args.concept_count_file,
         args.rs_entrez_symbol_file, args.prefix, args.no_norm_mutation_file)
