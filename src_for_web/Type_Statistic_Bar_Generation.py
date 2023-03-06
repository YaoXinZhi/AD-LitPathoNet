# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 31/05/2022 20:24
@Author: XINZHI YAO
"""

import os

import json

from collections import defaultdict

def entity_statistic(chain_file: str):

    type_to_id = defaultdict(set)  # 种类id统计

    reg_set = {'Reg', 'NegReg', 'PosReg'}

    prefix_to_type = {
        'GO': 'GO',
        'rs': 'Var',
        'HP': 'HPO',
        'OM': 'Disease',
        'ME': 'Disease'
    }
    with open(chain_file) as f:# 打开文件
        for line in f:
            l = line.strip().split('\t')  # 读取TXT数据，去掉多余空格并分割

            tags = l[3:]
            for tag in tags:
                if tag == 'ThemeOf' or tag == 'CauseOf':  # 消除对eval()干扰
                    continue

                # fixme: only normalized entity
                tag = eval(tag)  # 把数据存储为元组
                if len(tag) == 4 and tag[-1] != '-' and tag[1] not in reg_set:  # 去除Reg项
                    tag_type = tag[1]
                    tag_id = tag[-1].split('-')[-1]
                    if tag_type == 'Gene':  # 去掉'Gene:'
                        type_to_id[tag_type].add(tag_id)  # 将新id统计
                    else:

                        try:
                            nom_type = prefix_to_type[tag_id[:2]]  # 通过字典映射为种类
                        except:
                            # print(tag)
                            # print(tag_id)
                            # input()
                            continue
                        type_to_id[nom_type].add(tag_id)  # 统计映射后种类的id

    for _type, _ in type_to_id.items():
        'type_to_count[_type] = len(id_set)'

    type_to_count = {_type: len(id_set) for _type, id_set in
                     type_to_id.items()}  # 计数，等价于type_to_count[_type] = len(id_set)

    return type_to_count

def read_html(html_frag_file: str):
    html_doc = open(html_frag_file).read()
    return html_doc

def data_to_echarts(type_to_count: dict):
    data_str = """
    dataset: {
        dimensions: ['product', 'Count'],
        source: [ 
    """

    for _type in type_to_count.keys():
        if _type == 'Gene':
            counts = type_to_count[ _type ]
            node_str = "{" + \
                       'product: ' + "'Gene'" + " , 'Count': " + str(counts) + \
                       '},'
            data_str += node_str
        elif _type == 'Var':
            counts = type_to_count[ _type ]
            node_str = "{" + \
                       'product: ' + "'Var'" + " , 'Count': " + str(counts) + \
                       '},'
            data_str += node_str
        elif _type == 'GO':
            counts = type_to_count[ _type ]
            node_str = "{" + \
                       'product: ' + "'GO'" + " , 'Count': " + str(counts) + \
                       "},"
            data_str += node_str
        elif _type == 'HPO':
            counts = type_to_count[ _type ]
            node_str = "{" + \
                       'product: ' + "'HPO'" + " , 'Count': " + str(counts) + \
                       '},'
            data_str += node_str
        elif _type == 'Disease':
            counts = type_to_count[ _type ]
            node_str = "{" + \
                       'product: ' + "'Disease'" + " , 'Count': " + str(counts) + \
                       '},'
            data_str += node_str
    data_str += "]" + '\n' + '}'

    return data_str

def generate_type_sta_bar_option(chain_file: str,
                                 bar_html_file: str='../data/html_fragments/type_sta_bar_option.txt'):
    #dataset_replace_here
    bar_option = read_html(bar_html_file)

    type_to_count = entity_statistic(chain_file)

    data_str = data_to_echarts(type_to_count)

    bar_option = bar_option.replace('dataset_replace_here', data_str)

    return bar_option


# if __name__ == '__main__':
#
#     parser = argparse.ArgumentParser(description='Bar Echarts.')
#
#     parser.add_argument('-cf', dest='chain_file', required=True,
#                         help='chain file.')
#     parser.add_argument('-sf', dest='save_html', required=True,
#                         help='saved html file, *.html.')
#     args = parser.parse_args()
#
#     main(args.chain_file, args.save_html)


