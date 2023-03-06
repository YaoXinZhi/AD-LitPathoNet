# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 01/08/2022 14:45
@Author: yao
"""

import os
import re

import gzip
import argparse

from collections import defaultdict

def offset_adjust(entity1: str, sent: str, entity1_start: str, entity1_end: str):

    entity1_sent = sent[entity1_start: entity1_end]

    write_flag = False
    if entity1 == entity1_sent:
        write_flag = True
        return write_flag, entity1_start, entity1_end

    if entity1 == entity1_sent[ 1: ]:
        entity1_start += 1
        write_flag = True
    elif entity1 == entity1_sent[ 2: ]:
        entity1_start += 2
        write_flag = True
    elif entity1 == entity1_sent[ :-1 ]:
        entity1_end -= 1
        write_flag = True
    elif entity1 == entity1_sent[ :-2 ]:
        entity1_end -= 2
        write_flag = True
    elif entity1 == sent[ entity1_start - 1: entity1_end - 2 ]:
        entity1_start -= 1
        entity1_end -= 2
        write_flag = True
    elif entity1 == entity1_sent[ 1:-1 ]:
        write_flag = True
        entity1_start += 1
        entity1_end -= 1

    if entity1 != sent[entity1_start: entity1_end]:
        # wrong offset and can't adjust
        write_flag = False
    return write_flag, entity1_start, entity1_end

def offset_correction(src_file: str, Correction_file: str, reg_save_file=None):
    wf = open(Correction_file, 'w')

    wrong_count = 0
    total_count = 0
    reg_to_type = defaultdict(set)

    with open(src_file) as f:
        line = f.readline()
        wf.write(line)

        for line in f:
            l = line.strip().split('\t')
            entity1, entity1_type, entity1_offset, entity1_norm, entity1_agac = l[ :5 ]
            entity2, entity2_type, entity2_offset, entity2_norm, entity2_agac = l[ 7: 12 ]
            try:
                entity1_offset = eval(entity1_offset)
                entity2_offset = eval(entity2_offset)
            except:
                print(f'wrong offset parsing')
                print(entity1_offset)
                print(entity2_offset)
                continue

            entity1_start, entity1_end = map(int, entity1_offset)
            entity2_start, entity2_end = map(int, entity2_offset)

            entity1_end += 1
            entity2_end += 1

            sent = l[ -4 ]
            # pmid = l[ -3 ]
            # sent_id = l[ -2 ]
            # relation = l[ -1 ]

            # entity1_sent = sent[ entity1_start: entity1_end ]
            # entity2_sent = sent[ entity2_start: entity2_end ]
            total_count += 2

            # reg的标注保存
            if entity1_type in {'PosReg', 'NegReg', 'Reg'}:
                reg_to_type[ entity1 ].add(entity1_type)

            if entity2_type in {'PosReg', 'NegReg', 'Reg'}:
                reg_to_type[ entity2 ].add(entity2_type)

            entity1_flag, entity1_start, entity1_end = offset_adjust(entity1, sent, entity1_start, entity1_end)
            entity2_flag, entity2_start, entity2_end = offset_adjust(entity2, sent, entity2_start, entity2_end)

            entity1_offset = (entity1_start, entity1_end)
            entity2_offset = (entity2_start, entity2_end)

            l[ 2 ] = entity1_offset
            l[ 9 ] = entity2_offset

            if entity1_flag and entity2_flag:
                wf_line = '\t'.join(map(str, l))
                wf.write(f'{wf_line}\n')
            else:
                continue

        print(f'wrong offset: {wrong_count:,}/{total_count:,}')

    wf.close()
    print(f'{Correction_file} save done.')

    if not reg_save_file is None:
        with open(reg_save_file, 'w') as wf:
            wf.write(f'Reg\ttype\ttype_count\n')
            for reg, type_set in reg_to_type.items():
                type_wf = ', '.join(type_set)
                type_count = len(type_set)

                wf.write(f'{reg}\t{type_wf}\t{type_count}\n')
        print(f'{reg_save_file} save done.')


def read_homo_gene_file(homo_file: str):
    group_to_homo_id = {}

    gene_to_group = defaultdict()

    homo_entrez_to_symbol = {}
    with open(homo_file) as f:
        for line in f:
            l = line.strip().split('\t')

            group = l[ 0 ]
            tax_id = l[ 1 ]
            entrez = l[ 2 ]
            symbol = l[ 3 ]

            if tax_id == '9606':
                group_to_homo_id[ group ] = entrez
                homo_entrez_to_symbol[ entrez ] = symbol

            gene_to_group[ entrez ] = group

    filter_gene_to_group = {entrez: group for entrez, group in gene_to_group.items() if group_to_homo_id.get(group)}

    return filter_gene_to_group, group_to_homo_id, homo_entrez_to_symbol

# map entrez id to homo entrez
def entrez_mapping_to_homo(homo_entrez_to_symbol: dict, group_to_homo_id: dict, gene_to_group: dict, entrez: str):
    if homo_entrez_to_symbol.get(entrez):
        return entrez
    else:
        if gene_to_group.get(entrez):
            group = gene_to_group[ entrez ]

            homo_entrez = group_to_homo_id[ group ]
            return homo_entrez
        else:
            return ''

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

def read_umls(umls_file: str):
    id_to_name = {}

    with open(umls_file) as f:
        for line in f:
            l = line.strip().split('|')

            primary_log = l[ 2 ]

            _id = l[ 10 ]
            name = l[ 14 ]
            #             print(l)
            if primary_log == 'P':
                if id_to_name.get(_id):
                    if len(name) < len(id_to_name[ _id ]) and ',' not in name:
                        id_to_name[ _id ] = name
                else:
                    id_to_name[ _id ] = name
            else:
                if not id_to_name.get(_id):
                    id_to_name[ _id ] = name

    return id_to_name

def read_gene_info(gene_info_file: str):
    entrez_to_symbol = {}
    with gzip.open(gene_info_file) as f:
        f.readline()
        for line in f:
            l = line.decode('utf-8').strip().split('\t')
            entrez = l[ 1 ]
            symbol = l[ 2 ]

            entrez_to_symbol[ entrez ] = symbol

    print(f'gene count: {len(entrez_to_symbol):,}')
    return entrez_to_symbol

def get_new_entity_norm(entity_norm: str, homo_entrez_to_symbol: dict,
                        group_to_homo_id: dict, gene_to_group: dict):
    """
    转换为智人基因id
    """
    entity_new_tag_set = set()
    for tag in entity_norm.split(', '):
        tag_split = tag.split('-')

        if len(tag_split) < 3:
            continue

        if tag_split[ -2 ] == 'Gene':

            entrez = tag_split[ -1 ]
            if len(entrez.split(';')) > 1:
                entrez = entrez.split(';')[ 0 ]

            homo_entrez = entrez_mapping_to_homo(homo_entrez_to_symbol, group_to_homo_id, gene_to_group, entrez)
            if not homo_entrez:
                continue
            tag_split[ -1 ] = homo_entrez
            tag = '-'.join(tag_split)
            entity_new_tag_set.add(tag)
        else:
            entity_new_tag_set.add(tag)
    # print(entity_new_tag_set)

    if entity_new_tag_set:
        return ', '.join(entity_new_tag_set)
    else:
        return '-'


def get_nor_dis_name_from_umls(entity_norm: str, umls_id_to_name: dict):
    entity_new_tag_set = set()
    for tag in entity_norm.split(', '):
        tag_split = tag.split('-')

        if len(tag_split) < 3:
            continue

        if tag_split[ -2 ] == 'Disease':
            dis_id = tag_split[ -1 ]

            if dis_id.startswith('MESH'):

                _id = dis_id[ 5: ]

                if umls_id_to_name.get(_id):
                    umls_name = umls_id_to_name[ _id ]
                else:
                    continue

                tag_split[ 0 ] = umls_name
                tag = '-'.join(tag_split)
                entity_new_tag_set.add(tag)
            else:
                entity_new_tag_set.add(tag)
        else:
            entity_new_tag_set.add(tag)

    if entity_new_tag_set:
        return ', '.join(entity_new_tag_set)
    else:
        return '-'


def norm_to_norm_name(entity_norm: str, go_id_to_name: dict,
                      hpo_id_to_name: dict, entrez_to_symbol: dict,
                      go_to_type: dict):
    """
    1. 把GO HPO都换为标准名称
    2. 删除 CC的GO
    """
    new_norm_set = set()
    for norm in entity_norm.split(', '):
        norm_split = norm.split('-')
        norm_id = norm_split[ -1 ]
        norm_type = norm_split[ -2 ]

        if norm_id.startswith('GO:'):
            # print(norm)
            norm_name = go_id_to_name[ norm_id ]

            _type = go_to_type[ norm_id ] if go_to_type.get(norm_id) else ''

            if _type == 'cellular_component':
                continue

            if norm_name.startswith('obsolete'):
                continue
        elif norm_id.startswith('HP:'):
            if len(norm_id.split(';')) > 1:
                norm_id = norm_id.split(';')[ 0 ]
            norm_name = hpo_id_to_name[ norm_id ]
            if norm_name.startswith('obsolete'):
                continue
        elif norm_type == 'Gene':
            norm_name = entrez_to_symbol[ norm_id ]
        else:
            norm_name = ''
        if norm_name:
            new_norm_set.add(f'{norm_name}-{norm_type}-{norm_id}')
        else:
            new_norm_set.add(norm)
    if new_norm_set:
        return ', '.join(new_norm_set)
    else:
        return '-'

def agac_label_convert(label: str):
    label_convert = {
        # convert to Gene
        'Protein': 'Gene',
        'Enzyme': 'Gene',
        # convert to CPA/MPA
        'Interaction': 'CPA',
        'Pathway': 'CPA'
    }
    if label_convert.get(label):
        return label_convert[label]
    else:
        return label

def delete_wrong_type_norm(entity_norm: str, entity_type: str):
    process_type_set = {'Disease', 'Pathway', 'CPA', 'MPA', 'Interaction'}

    new_norm_set = set()
    for norm in entity_norm.split(', '):

        if norm == '-' or len(norm.split('-')) < 2:
            continue

        norm_split = norm.split('-')
        # print(norm_split)
        norm_id = norm_split[ -1 ]
        norm_type = norm_split[ -2 ]

        if norm_type in {'Chemical', 'Species'}:
            continue

        if entity_type == 'Gene':
            if norm_type not in {'Gene', 'cellular_component'}:
                continue

        if entity_type == 'Var':
            if norm_type != 'Mutation':
                continue

        if entity_type in {'Protein', 'Enzyme'}:
            if norm_type not in {"cellular_component", 'Gene'}:
                continue

        if entity_type in process_type_set:
            if norm_id[ :2 ] not in {'GO', 'HP', 'ME', 'OM'}:
                continue

        new_norm_set.add(norm)

    if not new_norm_set:
        return '-'
    else:
        return ', '.join(new_norm_set)

def only_max_length(entity_norm: str):
    norm_list = entity_norm.split(', ')

    if len(norm_list) <= 1:
        return entity_norm

    max_len_norm = [ norm.split('-')[ 0 ] for norm in norm_list ]
    max_len_norm.sort()
    max_len_norm = max_len_norm[ -1 ]
    new_norm_set = set()
    for norm in norm_list:
        norm_split = norm.split('-')
        norm_mention = norm_split[ 0 ]

        if norm_mention == max_len_norm:
            new_norm_set.add(norm)

    if not new_norm_set:
        print('Error!!!')
        print(new_norm_set)
        input()

    return ', '.join(new_norm_set)

def gene_to_homo(src_file: str, Correction_file: str,
                 homo_entrez_to_symbol: dict, group_to_homo_id: dict,
                 gene_to_group: dict, homo_symbol_to_entrez: dict,
                 go_id_to_name: dict, hpo_id_to_name: dict, entrez_to_symbol: dict,
                 umls_id_to_name,
                 #no_annotation_go_set: set,
                 go_to_type: dict,
                 missed_save_file: str = ''):
    """
    1. 将其他物种的entrez id 对应到人类的entrez id
    2. 将缩写明确的 直接对应到人类entrez id
    """

    wf = open(Correction_file, 'w')
    missed_set = set()
    with open(src_file) as f:
        line = f.readline()
        wf.write(line)

        for line in f:
            l = line.strip().split('\t')
            entity1, entity1_type, entity1_offset, entity1_tool, entity1_agac, entity1_norm, entity1_norm_count = l[:7]
            entity2, entity2_type, entity2_offset, entity2_tool, entity2_agac, entity2_norm, entity2_norm_count = l[7: 14 ]

            # Convert agac label for higher recall
            entity1_type = agac_label_convert(entity1_type)
            entity2_type = agac_label_convert(entity2_type)

            l[1] = entity1_type
            l[8] = entity2_type

            # 匹配同源的人类基因
            new_entity1_norm = get_new_entity_norm(entity1_norm, homo_entrez_to_symbol,
                                                   group_to_homo_id, gene_to_group)

            new_entity2_norm = get_new_entity_norm(entity2_norm, homo_entrez_to_symbol,
                                                   group_to_homo_id, gene_to_group)

            # 给疾病匹配标准的名称 目前只安排的umls的
            new_entity1_norm = get_nor_dis_name_from_umls(new_entity1_norm, umls_id_to_name)

            new_entity2_norm = get_nor_dis_name_from_umls(new_entity2_norm, umls_id_to_name)

            # 匹配人类基因的symbol
            if entity1_type == 'Gene' and (entity1_norm == '-' or new_entity1_norm == '-'):
                if homo_symbol_to_entrez.get(entity1.upper()):
                    homo_entrez = homo_symbol_to_entrez[ entity1.upper() ]
                    new_entity1_norm = f'{entity1}-Gene-{homo_entrez}'
                # else:
                # missed_set.add((entity1, entity1_norm))
            if entity2_type == 'Gene' and (entity2_norm == '-' or new_entity2_norm == '-'):
                if homo_symbol_to_entrez.get(entity2.upper()):
                    homo_entrez = homo_symbol_to_entrez[ entity2.upper() ]
                    new_entity2_norm = f'{entity2}-Gene-{homo_entrez}'
                # else:
                # missed_set.add((entity2, entity2_norm))

            # 给rs开头的突变加上标准化
            if entity1_type == 'Var' and (entity1_norm == '-' or new_entity1_norm == '-'):
                if entity1.lower().startswith('rs'):
                    rs_id = re.findall(r'rs[\d]+', entity1.lower())
                    if rs_id:
                        rs_id = rs_id[ 0 ]
                        new_entity1_norm = f'{rs_id}-Mutation-{rs_id}'

            if entity2_type == 'Var' and (entity2_norm == '-' or new_entity2_norm == '-'):
                if entity2.lower().startswith('rs'):
                    rs_id = re.findall(r'rs[\d]+', entity2.lower())
                    if rs_id:
                        rs_id = rs_id[ 0 ]
                        new_entity2_norm = f'{rs_id}-Mutation-{rs_id}'

            # 删除和实体类型不一致的标准化id
            new_entity1_norm = delete_wrong_type_norm(new_entity1_norm, entity1_type)
            new_entity2_norm = delete_wrong_type_norm(new_entity2_norm, entity2_type)

            # 仅仅保留最长的entity_norm
            # 在删除错误的标注类型之后
            new_entity1_norm = only_max_length(new_entity1_norm)
            new_entity2_norm = only_max_length(new_entity2_norm)

            #  给norm列的标注 entrez，GO，HPO都换为标准name
            new_entity1_norm = norm_to_norm_name(new_entity1_norm, go_id_to_name,
                                                 hpo_id_to_name, entrez_to_symbol,
                                                 go_to_type)
            new_entity2_norm = norm_to_norm_name(new_entity2_norm, go_id_to_name,
                                                 hpo_id_to_name, entrez_to_symbol,
                                                 go_to_type)

            l[ 5 ] = new_entity1_norm
            l[ 12 ] = new_entity2_norm

            wf_line = '\t'.join(l)

            wf.write(f'{wf_line}\n')

    wf.close()
    print(f'{Correction_file} save done.')

    if missed_save_file:
        with open(missed_save_file, 'w') as wf:
            for mention, norm in missed_set:
                wf.write(f'{mention}\t{norm}\n')

    print(f'{missed_save_file} save done.')

def main(infer_merge_file: str, save_path: str, prefix: str):

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    # processing file
    correction_file = f'{save_path}/{prefix}.correction.tsv'
    reg_save_file = f'{save_path}/{prefix}.reg.tsv'

    final_file = f'{save_path}/{prefix}.correction-homo.tsv'
    missed_save_file = f'{save_path}/{prefix}.miss.tsv'

    # some source file.
    homo_gene_mapping_file = '../data/ncbi-data/homologene.data.txt'
    mesh_id_to_str_file = '../data/umls/MRCONSO.RRF'

    go_file = '../data/obo-data/go.obo'

    hpo_file = '../data/obo-data/hp.obo'

    homo_info_file = '../data/ncbi-data/Homo_sapiens.gene_info.gz'

    print(f'1. Offset correction.')
    offset_correction(infer_merge_file, correction_file, reg_save_file)

    print(f'2. Reading database file.')
    gene_to_group, group_to_homo_id, homo_entrez_to_symbol = read_homo_gene_file(homo_gene_mapping_file)

    umls_id_to_name = read_umls(mesh_id_to_str_file)
    go_id_to_name, no_annotation_go_set, go_to_type = read_obo(go_file)

    hpo_id_to_name, _, _ = read_obo(hpo_file)

    entrez_to_symbol_ncbi = read_gene_info(homo_info_file)

    # 合并两个来源的entrez_to_symbol (homo_mapping, ncbi_gene_info)
    entrez_to_symbol = {}

    for entrez, symbol in entrez_to_symbol_ncbi.items():
        entrez_to_symbol[entrez] = symbol

    for entrez, symbol in homo_entrez_to_symbol.items():
        if not entrez_to_symbol.get(entrez):
            entrez_to_symbol[entrez] = symbol
    #print(f'total gene count: {len(entrez_to_symbol):,}')

    homo_symbol_to_entrez = {symbol: entrez for entrez, symbol in entrez_to_symbol.items()}

    print(f'3. merge Homo gene info.')
    gene_to_homo(src_file=correction_file, Correction_file=final_file,
                 homo_entrez_to_symbol=homo_entrez_to_symbol,
                 group_to_homo_id=group_to_homo_id,
                 gene_to_group=gene_to_group,
                 homo_symbol_to_entrez=homo_symbol_to_entrez,
                 go_id_to_name=go_id_to_name, hpo_id_to_name=hpo_id_to_name,
                 entrez_to_symbol=entrez_to_symbol, umls_id_to_name=umls_id_to_name,
                 go_to_type=go_to_type, missed_save_file=missed_save_file
                 )





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Entity Normaliztion.')
    parser.add_argument('-if', dest='infer_merge_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    parser.add_argument('-pr', dest='prefix', required=True)
    args = parser.parse_args()

    main(args.infer_merge_file, args.save_path, args.prefix)
