# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 11/11/2022 21:03
@Author: XINZHI YAO
"""

"""
改代码用于基础版本的数据库合并
因为大数据版本数据居然更少了？？
"""

import os

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


def read_file(_file: str, mutation_type_set: set, mutation_sub_type_set: set):
    print(f'Read: {_file}')

    mutation_type_idx = ''
    sub_mutation_idx = ''

    l_set = set()
    with open(_file) as f:

        head = f.readline().strip()

        print(mutation_type_set)
        if not _file.endswith('evi.tsv'):
            mutation_type_idx = head.split('\t').index('MutationType')
            sub_mutation_idx = head.split('\t').index('MutationSubType')
            chain_index = head.split('\t').index('Chain_line')

        for line in f:
            l = line.strip()

            if mutation_type_idx:
                l_split = l.split('\t')
                mutation_type = l_split[mutation_type_idx]
                mutation_sub_type = l_split[sub_mutation_idx]

                if (mutation_type not in mutation_type_set or mutation_sub_type not in mutation_sub_type_set) \
                        and mutation_sub_type!= 'RsMutation':

                    chain = l_split[chain_index:]

                    l_split[mutation_type_idx] = 'point mutations'
                    l_split[sub_mutation_idx] = 'RsMutation'

                    for idx, ele in enumerate(chain):
                        if ", 'Var', " in ele:
                            var_ele = ele
                            var_idx = chain_index + idx
                            break
                    try:
                        var_ele = eval(var_ele)

                        var_ele = (var_ele[0], var_ele[1], var_ele[2], 'point mutations-Mutation-Mutations')

                        l_split[var_idx] = str(var_ele)

                        l = '\t'.join(l_split)
                    except:
                        print(mutation_type in mutation_type_set)
                        print(mutation_sub_type in mutation_sub_type_set)
                        print(mutation_type)
                        print(mutation_sub_type)
                        print(f'var_ele {var_ele}')
                        print(f'var_idx_ele: {l_split[var_idx]}')
                        print(l_split[chain_index:])
                        continue
            l_set.add(l)
    return head, l_set


def merge_db_file(file_a: str, file_b: str, save_file: str, mutation_type_set: set, sub_mutation_type_set: set):

    save_set = set()

    head, l_set = read_file(file_a, mutation_type_set, sub_mutation_type_set)
    save_set.update(l_set)

    _, _l_set = read_file(file_b, mutation_type_set, sub_mutation_type_set)
    save_set.update(_l_set)

    with open(save_file, 'w') as wf:
        wf.write(f'{head}\n')
        for l in save_set:
            wf.write(f'{l}\n')

    print(f'{save_file} save done, {len(save_set):,} line saved.')

base_path = '/home/xzyao/AD-PNRLE/result/complete_database'

no_norm_mutation = '/home/xzyao/AD-PNRLE/data/Non-Normalized_Mutations/MutationTypeDictionary20220711.txt'

big_prefix = 'total.filter'
big_db_path = f'{base_path}/base_DB_dir'

small_prefix = 'ad.small'
small_db_path = f'{base_path}/small_base_dir'

save_prefix = 'iASiS-PNRLE'
save_path = f'{base_path}/base_db_merge_dir'

if not os.path.exists(save_path):
    os.mkdir(save_path)

merge_type_list = ['DjangoDB', 'GO', 'HPO', 'PPI', 'evi', 'DIS']

mutation_set, sub_mutation_set = read_no_norm_mutation_file(no_norm_mutation)

for merge_type in merge_type_list:
    print(f'Processing: {merge_type}.')
    small_db = f'{small_db_path}/{small_prefix}.{merge_type}.tsv'
    big_db = f'{big_db_path}/{big_prefix}.{merge_type}.tsv'

    save_db = f'{save_path}/{save_prefix}.{merge_type}.tsv'

    merge_db_file(small_db, big_db, save_db, mutation_set, sub_mutation_set)
