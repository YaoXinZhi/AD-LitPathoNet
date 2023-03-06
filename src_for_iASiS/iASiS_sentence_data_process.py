# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 05/12/2022 10:14
@Author: XINZHI YAO
"""

"""
该代码用于筛选哪些iASiS句子需要存到数据库
只讲iASiS作为额外的边在已有的PNRLE上进行添加
1. 统计已有的所有表型条目
2. 筛选所有的iASiS的边 （Functional Change）
3. 晒徐娜所有关联到的句子 PMID 存成数据库格式
"""

import sys
import json
import re
import os
from collections import defaultdict

def all_concept_stat(concept_db_file: str):

    concept_id_set = set()
    with open(concept_db_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            term_id = l[2]

            if term_id.startswith('MESH'):
                term_id = term_id.split(':')[1]

            # term_name = l[3]

            concept_id_set.add(term_id)

    print(f'total {len(concept_id_set):,} concepts.')
    return concept_id_set

def iasis_sent_id_stat(concept_set: set, iasis_funrel_file: str):

    used_sent_id_set = set()
    match_concept_pair_count = 0

    triple_to_sent_id = defaultdict(set)
    with open(iasis_funrel_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            relation = l[0]

            concept_id_1 = l[3]
            concept_id_2 = l[7]

            sent_id = re.findall(r"sent_id: \[(.*?)]", l[ 1 ])

            # print(sent_id)
            if len(sent_id) < 1:
                continue

            sent_id = sent_id[0]

            sent_id_list = re.findall(r"'(.*?)'", sent_id)

            # if len(sent_id_list) > 1:
            #     print(sent_id_list)
            #     input()

            if concept_id_1 in concept_set and concept_id_2 in concept_set:
                match_concept_pair_count += 1

                used_sent_id_set.update(sent_id_list)

                triple_to_sent_id[(concept_id_1, relation, concept_id_2)].update(sent_id_list)

    print(f'match_concept_pair: {match_concept_pair_count:,}')
    print(f'used_sent_count: {len(used_sent_id_set):,}')
    return used_sent_id_set, triple_to_sent_id

def json_process(json_doc: json, mode: str):

    pmid = json_doc['id']
    if 'NumberInt' in pmid:
        pmid = re.findall(r'NumberInt\((\d+)\)', pmid)[0]

    sent_id_to_text = {}
    sent_list = json_doc['sents']
    for sent_json in sent_list:
        sent_id = sent_json['sent_id']

        sent_id = f'{pmid}_{mode}_{sent_id}'

        sent_text = sent_json['sent_text']

        sent_id_to_text[sent_id] = sent_text

    return sent_id_to_text

def iasis_sent_file_process_report_err(iasis_json_file: str, sent_id_set: set):

    if 'Abstracts' in iasis_json_file:
        mode = 'abstract'
    else:
        mode = 'fullText'

    doc = ''
    save_sent_id_to_text = {}
    process_count = 0
    with open(iasis_json_file) as f:
        for line in f:
            process_count += 1
            if process_count%500000 == 0:
                print(f'{process_count:,} processed.')

            if line[ 0 ] == '{':
                if not doc:
                    doc += line
                    continue
                else:
                    doc = ''
                    doc += line

            elif line[ 0 ] == '}':
                doc += line
                json_doc = eval(doc)
                # json_list.append(json_doc)

                sent_id_to_text = json_process(json_doc, mode)

                for sent_id, text in sent_id_to_text.items():
                    if sent_id in sent_id_set:
                        save_sent_id_to_text[sent_id] = text

            else:
                if 'NumberInt' in line:
                    # print(line)
                    try:
                        id_str = re.findall(r'NumberInt\(\d+\)', line)[ 0 ]
                    except:
                        id_str = re.findall(r'NumberInt\(.*?\)', line)[ 0 ]
                    line = line.replace(id_str, f'"{id_str}"')

                if 'ObjectId' in line:
                    # print(line)
                    id_str = re.findall(r'ObjectId\(.*?\)', line)[ 0 ]
                    new_id_str = id_str.replace('"', "'")
                    line = line.replace(id_str, f'"{new_id_str}"')

                doc += line

    print(f'{len(save_sent_id_to_text):,} sentence saved.')
    return save_sent_id_to_text

def iasis_sent_file_process(iasis_json_file: str, sent_id_set: set):

    if 'Abstracts' in iasis_json_file:
        mode = 'abstract'
    else:
        mode = 'fullText'

    doc = ''
    save_sent_id_to_text = {}
    process_count = 0
    try:
        with open(iasis_json_file) as f:
            for line in f:
                process_count += 1
                if process_count%500000 == 0:
                    print(f'{process_count:,} processed.')

                if line[ 0 ] == '{':
                    if not doc:
                        doc += line
                        continue
                    else:
                        doc = ''
                        doc += line

                elif line[ 0 ] == '}':
                    doc += line
                    json_doc = eval(doc)
                    # json_list.append(json_doc)

                    sent_id_to_text = json_process(json_doc, mode)

                    for sent_id, text in sent_id_to_text.items():
                        if sent_id in sent_id_set:
                            save_sent_id_to_text[sent_id] = text

                else:
                    if 'NumberInt' in line:
                        # print(line)
                        id_str = re.findall(r'NumberInt\(\d+\)', line)[ 0 ]
                        line = line.replace(id_str, f'"{id_str}"')

                    if 'ObjectId' in line:
                        # print(line)
                        id_str = re.findall(r'ObjectId\(.*?\)', line)[ 0 ]
                        new_id_str = id_str.replace('"', "'")
                        line = line.replace(id_str, f'"{new_id_str}"')

                    doc += line
    except:
        print(sys.exc_info())
        input()
        with open(iasis_json_file, encoding='unicode_escape') as f:
            for line in f:
                process_count += 1
                if process_count % 500000 == 0:
                    print(f'{process_count:,} processed.')

                if line[ 0 ] == '{':
                    if not doc:
                        doc += line
                        continue
                    else:
                        doc = ''
                        doc += line

                elif line[ 0 ] == '}':
                    doc += line
                    try:
                        json_doc = eval(doc)
                    except:
                        print(sys.exc_info())
                        print(doc)
                        input()
                    # json_list.append(json_doc)

                    sent_id_to_text = json_process(json_doc, mode)

                    for sent_id, text in sent_id_to_text.items():
                        if sent_id in sent_id_set:
                            save_sent_id_to_text[ sent_id ] = text

                else:
                    if 'NumberInt' in line:
                        # print(line)
                        id_str = re.findall(r'NumberInt\(\d+\)', line)[ 0 ]
                        line = line.replace(id_str, f'"{id_str}"')

                    if 'ObjectId' in line:
                        # print(line)
                        id_str = re.findall(r'ObjectId\(.*?\)', line)[ 0 ]
                        new_id_str = id_str.replace('"', "'")
                        line = line.replace(id_str, f'"{new_id_str}"')

                    doc += line

    print(f'{len(save_sent_id_to_text):,} sentence saved.')
    return save_sent_id_to_text

def batch_iasis_json(iasis_json_file_list: list, sent_id_set: set):
    total_save_sent_id_to_text = {}


    for _file in iasis_json_file_list:
        print(f'Processing: {_file}')
        save_sent_id_to_text = iasis_sent_file_process(_file, sent_id_set)
        # save_sent_id_to_text = iasis_sent_file_process_report_err(_file, sent_id_set)

        total_save_sent_id_to_text.update(save_sent_id_to_text)
    return total_save_sent_id_to_text

def save_iasis_evi_db_file(triple_to_sent_id: dict,
                           total_save_sent_id_to_text: dict, iasis_evi_save_file: str):

    match_count = 0
    miss_set = set()
    save_count = 0
    with open(iasis_evi_save_file, 'w') as wf:
        wf.write('Relation\tTerm1\tTerm2\tPMID\tSentID\tSentence\n')
        for triple, sent_id_set in triple_to_sent_id.items():
            term_1, relation, term_2 = triple
            for sent_id in sent_id_set:
                pmid = sent_id.split('_')[0]
                if total_save_sent_id_to_text.get(sent_id):
                    sentence = total_save_sent_id_to_text[sent_id ]
                    match_count += 1
                else:
                    miss_set.add(sent_id)
                    # sentence = 'None'
                    continue

                wf.write(f'{relation}\t{term_1}\t{term_2}\t'
                         f'{pmid}\t{sent_id}\t'
                         f'{sentence.strip()}\n')
            save_count += 1

    print(miss_set)
    print(f'{match_count:,} sent_id matched.')
    print(f'{len(miss_set):,} sent_id missed.')
    print(f'{iasis_evi_save_file} save done, {save_count:,} evidence saved.')


def main():

    Bioconcept_db_file = '/mnt/disk1/xzyao/AD-PNRLE/result/base_db_merge_dir/PNRLE.BioConcept.tsv'

    # iASiS concept id, relation type, sentence id et.al.
    iASiS_funrel_file = '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_all_convert/PNRLE_from_iASiS.all.tsv'
    # iASiS_funrel_file = '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_all_convert/PNRLE_from_iASiS.all.FunRel.tsv'

    iasis_json_file_list = [
        '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_Sentence_data/sentences_a/AD_Abstracts_ENRICHED.json',
        '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_Sentence_data/sentences_a/AD_FullText_ENRICHED.json',

        '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_Sentence_data/sentences_b/AD_Abstract_enriched.json',
        '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_Sentence_data/sentences_b/AD_Fulltext_enriched.json',

        # '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_Sentence_data/sentences_c/AD_Abstract_enriched.json',
        # '/mnt/disk1/xzyao/AD-PNRLE/data/iASiS_Sentence_data/sentences_c/AD_Fulltext_enriched.json'
    ]

    # iasis_evi_save_file = f'/mnt/disk1/xzyao/AD-PNRLE/result/base_db_merge_dir/iASiS.Evi.db.tsv'
    iasis_evi_save_file = f'/mnt/disk1/xzyao/AD-PNRLE/result/base_db_merge_dir/iASiS.Evi.db.all.tsv'

    concept_id_set = all_concept_stat(Bioconcept_db_file)

    sent_id_set, triple_to_sent_id = iasis_sent_id_stat(concept_id_set, iASiS_funrel_file)

    total_save_sent_id_to_text = batch_iasis_json(iasis_json_file_list, sent_id_set)

    save_iasis_evi_db_file(triple_to_sent_id, total_save_sent_id_to_text, iasis_evi_save_file)








if __name__ == '__main__':
    pass
    # main()


