# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 17/05/2022 9:53
@Author: XINZHI YAO
"""
import itertools
import os
import argparse

from collections import defaultdict
from itertools import product

"""
该代码通过阅读SemanticRoleChain文件生成Echarts-Graph的HTML
1. 不包含筛选功能。
2. 会生成临时中间文件对SemanticRoleChain进行简单处理。
3. SemanticRoleChain文件包含 PMID, Sent_id, Sent, Chain
"""

def chain_process(chain: list):

    reg_set = {'Reg', 'NegReg', 'PosReg'}

    target_set = {'MPA', 'CPA', 'Disease'}

    edge_set = set()
    cause_idx = chain.index('CauseOf')

    node_to_norm = defaultdict(set)
    causer = ''
    tagger = ''
    norm_id_to_name = {}

    i = 0
    while True:
        if i+3 < len(chain):
            triple = chain[i: i+3]
        else:
            triple = chain[i: len(chain)]

        if len(triple) != 3:
            print(chain)
            print(i)
            print(triple)
            raise ValueError
        source, relation, target = triple
        # print(source)
        # print(relation)
        # print(target)
        source = eval(source)
        target = eval(target)

        source_mention, source_type, _, source_norm = source
        target_mention, target_type, _, target_norm = target

        # delete chains containing unnormalized node
        if (source_norm == '-' and source_type not in reg_set) \
            or (target_norm == '-' and target_type not in reg_set):
            return set(), {}, {}

        # delete chains containing Alzheimer's Disease (MESH:D000544, HP:0002511)
        if (('MESH:D000544' in source_norm or 'HP:0002511' in source_norm) and source_type not in reg_set) \
                or (('MESH:D000544' in target_norm or 'HP:0002511' in target_norm) and target_type not in reg_set):
            return set(), {}, {}

        if relation == 'ThemeOf':
            # ThemeOf on the left of CauseOf
            if i < cause_idx:
                edge_set.add((source, ('ThemeOf', 'ThemeOf'), target))

                # Gene ThemeOf Var
                if source_type == 'Gene' and target_type == 'Var':
                    if len(source_norm.split(', '))> 1:
                        # print(source_norm)
                        # input()
                        # raise ValueError(f'source_norm: {source_norm}')
                        source_norm = source_norm.split(', ')[0]
                    entrez = source_norm.split('-')[-1]
                    symbol = source_norm.split('-')[0]

                    node_to_norm[source].add(f'Gene:{entrez}')
                    norm_id_to_name[f'Gene:{entrez}'] = symbol
                    if node_to_norm.get(target):
                        node_to_norm[target] = set()
                    target_id = target_norm.split('-')[-1]
                    if target_mention.startswith('rs'):
                        node_to_norm[target].add(f'{entrez}-Var:{target_mention}')
                        norm_id_to_name[f'{entrez}-Var:{target_mention}'] = f'{symbol}:{target_mention}'
                    elif target_id.startswith('rs'):
                        node_to_norm[target].add(f'{entrez}-Var:{target_id}')
                        norm_id_to_name[f'{entrez}-Var:{target_id}'] = f'{symbol}:{target_id}'
                    # no normalization
                    else:
                        pass

                # Gene ThemeOf GO/HPO/Disease et.al.
                if source_type == 'Gene' and target_type in target_set:
                    edge_set.add((source, ('ThemeOf', 'ThemeOf'), target))
                    entrez = source_norm.split('-')[-1]
                    symbol = source_norm.split('-')[0]

                    node_to_norm[ source ].add(f'Gene:{entrez}')
                    norm_id_to_name[f'Gene:{entrez}'] = symbol

                    for norm in target_norm.split(', '):
                        norm_id = norm.split('-')[-1]
                        norm_name = '-'.join(norm.split('-')[0:-2]).strip().capitalize()

                        if not (norm_id.startswith('GO:')
                                or norm_id.startswith('HP:')
                                # 6/6 delete OMIM tagging temporarily
                                # or norm_id.startswith('OMIM:')
                                or norm_id.startswith("MESH:")):
                            continue

                        if node_to_norm.get(target):
                            node_to_norm[target] = set()
                        node_to_norm[target].add(f'{entrez}-{norm_id}')
                        norm_id_to_name[f'{entrez}-{norm_id}'] = f'{symbol}:{norm_name}'

            # ThemeOf on the right of CauseOf
            else:
                # * CauseOf GO/HPO/Disease et.al.
                if source == tagger:
                    edge_set.add((causer, ((tagger[0], tagger[2]), tagger[1]), target))
                    for norm in target_norm.split(','):
                        norm_id = norm.split('-')[-1]
                        norm_name = '-'.join(norm.split('-')[:-2]).strip().capitalize()
                        if not (norm_id.startswith('GO:')
                                or norm_id.startswith('HP:')
                                # 6/6 delete OMIM tagging temporarily
                                # or norm_id.startswith('OMIM:')
                                or norm_id.startswith('MESH:')):
                            continue
                        node_to_norm[target].add(f'{norm_id}')
                        norm_id_to_name[f'{norm_id}'] = norm_name
                else:
                    # Gene ThemeOf GO/HPO/Disease et.al.
                    if target_type == 'Gene' and source_type in target_set:
                        edge_set.add((target, ('ThemeOf', 'ThemeOf'), source))
                        entrez = target_norm.split('-')[-1]
                        symbol = target_norm.split('-')[ 0 ]

                        node_to_norm[target].add(f'Gene:{entrez}')
                        norm_id_to_name[ f'Gene:{entrez}' ] = symbol

                        for norm in source_norm.split(', '):
                            norm_id = norm.split('-')[-1]
                            norm_name = '-'.join(norm.split('-')[:-2]).strip().capitalize()
                            if not (norm_id.startswith('GO:')
                                    or norm_id.startswith('HP:')
                                    # 6/6 delete OMIM tagging temporarily
                                    # or norm_id.startswith('OMIM:')
                                    or norm_id.startswith('MESH:')):
                                continue

                            if node_to_norm.get(source):
                                node_to_norm[source] = set()
                            node_to_norm[source].add(f'{entrez}-{norm_id}')
                            norm_id_to_name[f'{entrez}-{norm_id}'] = f'{symbol}:{norm_name}'
        # * CauseOf *
        elif relation == 'CauseOf':
            causer = source
            tagger = target
            # print(causer)
            # print(tagger)
            # input()
            # Var with no normalization
            if not node_to_norm[source]:

                if source_type == 'Var':
                    source_id = source_norm.split('-')[ -1 ]
                    if source_mention.startswith('rs'):
                        node_to_norm[ source ].add(f'Var:{source_mention}')
                        norm_id_to_name[f'Var:{source_mention}'] = f'{source_mention}'
                    elif source_id.startswith('rs'):
                        node_to_norm[ source ].add(f'Var:{source_id}')
                        norm_id_to_name[f'Var:{source_id}'] = f'{source_id}'
                    else:
                        # no normalization
                        pass
        else:
            raise ValueError(f'relation: {relation}')

        i+= 2
        if i >= len(chain)-1:
            break
    # print(edge_set)
    # input()

    return edge_set, node_to_norm, norm_id_to_name

def edge_to_norm(edge_set: set, node_to_norm: dict, evidence: set):
    # normalized node and edge.

    norm_edge_set = set()
    norm_node_set = set()


    norm_edge_to_evidence = {}
    pmid, sent_id, sent = evidence
    for edge in edge_set:
        source, (tagger_word, relation), target = edge
        # print(source, tagger_word, relation, target)
        # input()

        source_norm_set = node_to_norm[source]
        target_norm_set = node_to_norm[target]

        norm_node_set.update(node_to_norm[source])
        norm_node_set.update(node_to_norm[target])

        for norm_source, norm_target in itertools.product(source_norm_set, target_norm_set):
            norm_edge = (norm_source, relation, norm_target)
            norm_edge_set.add(norm_edge)
            norm_edge_to_evidence[norm_edge] = [source, target, (pmid, sent_id, sent, tagger_word)]

    return norm_edge_set, norm_node_set, norm_edge_to_evidence

def save_temp_node_edge_file(save_file: str, count_dic: dict):

    element_sort = sorted(count_dic, key=lambda x: count_dic[ x ], reverse=True)

    with open(save_file, 'w') as wf:
        for element in element_sort:
            count = count_dic[element]
            wf.write(f'{element}\t{count}\n')
    print(f'Temporary file: {save_file} save done.')

def evidence_select(edge_to_evidence: dict):
    """
    该代码用于证据筛选
    1. 首先筛选非ThemeOf边的证据，（Reg/PosReg/NegReg）
    2. 给ThemeOf的边提供证据，尽量确保提供的是在调控边里出现的文本证据
    """

    edge_to_selected_evidence = {}
    selected_id_set = set()
    # 1. select evidence for Reg/PosReg/NegReg edge.
    for edge, evidence_set in edge_to_evidence.items():
        # print(edge)
        # print(list(evidence_set)[1])
        # input()
        relation = edge[1]
        if relation == 'ThemeOf':
            continue

        # fixme: select the shortest evidence.
        selected_evidence = list(sorted(list(evidence_set), key=lambda x: len(x[2][2])))[0]

        pmid, sent_id, sentence = selected_evidence

        edge_to_selected_evidence[edge] = selected_evidence
        selected_id_set.add((pmid, sent_id))

    # 2. select evidence for ThemeOf edge.
    for edge, evidence_set in edge_to_evidence.items():
        relation = edge[1]

        if relation != 'ThemeOf':
            continue
        theme_evi_saved = False
        for evidence in evidence_set:
            source, target, (pmid, sent_id, sentence, tagger_word) = evidence
            if (pmid, sent_id) in selected_id_set:
                theme_evi_saved = True
                edge_to_selected_evidence[edge] = evidence

        if not theme_evi_saved:
            edge_to_selected_evidence[edge] = list(sorted(list(evidence_set), key=lambda x: len(x[2][2])))[0]

    # print(edge_to_selected_evidence)

    return edge_to_selected_evidence


def save_tmp_rich_file(edge_to_rich: dict, save_file: str):

    with open(save_file, 'w') as wf:
        wf.write(f'Edge\tRiceText\n')
        for edge, rich_text in edge_to_rich.items():
            wf.write(f'{edge}\t{rich_text}\n')


def read_chain_file(chain_file: str, temp_dir: str='../temp_dir'):

    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)

    edge_count = defaultdict(int)
    node_count = defaultdict(int)
    total_norm_id_to_name = {}
    total_norm_edge_to_evidence = defaultdict(list)
    with open(chain_file) as f:
        for line in f:
            l = line.strip().split('\t')
            # print(l)
            pmid, sent_id, sent = l[0], l[1], l[2]
            chain = l[3:]
            edge_set, node_to_norm, norm_id_to_name = chain_process(chain)

            norm_edge_set, norm_node_set, norm_edge_to_evidence = edge_to_norm(edge_set, node_to_norm, (pmid, sent_id, sent))
            # fixme norm_edge_to_evidence.value()
            #  [source, target, (pmid, sent_id, sent, tagger_word)]

            total_norm_id_to_name.update(norm_id_to_name)

            for norm_edge, evidence in norm_edge_to_evidence.items():

                total_norm_edge_to_evidence[norm_edge].append(evidence)

            # edge statistics
            for norm_edge in norm_edge_set:
                edge_count[norm_edge] += 1

            for node in norm_node_set:
                node_count[node] += 1

    edge_to_selected_evidence = evidence_select(total_norm_edge_to_evidence)
    # print(edge_to_selected_evidence)
    # input()
    edge_to_rich_evidence = get_html_rich(edge_to_selected_evidence)

    node_save_file = f'{temp_dir}/node.tsv'
    edge_save_file = f'{temp_dir}/edge.tsv'
    save_temp_node_edge_file(node_save_file, node_count)
    save_temp_node_edge_file(edge_save_file, edge_count)

    rich_evidence_file = f'{temp_dir}/edge.rich.tsv'
    save_tmp_rich_file(edge_to_rich_evidence, rich_evidence_file)

    return node_count, edge_count, total_norm_id_to_name, edge_to_rich_evidence

def id_to_type(node_id: str):
    if 'Gene' in node_id:
        node_type = 'Gene'
    elif 'Var' in node_id:
        node_type = 'Var'
    elif 'GO' in node_id:
        node_type = 'GO'
    elif 'HP' in node_id:
        node_type = 'HPO'
    elif 'OMIM' in node_id:
        node_type = 'Disease'
    elif 'MESH' in node_id:
        node_type = 'Disease'
    else:
        node_type = ''
    return node_type

def get_html_rich(edge_to_selected_evidence: dict):

    type_to_color = {
        'Gene': '#ffa235',
        'Var': '#8932a5',
        'GO': '#90ed7d',
        'HPO': '#70ad47',
        'Disease': '#82B0D2',

        'Reg': '#DCDCDC',
        'PosReg': '#FF6666',
        'NegReg': '#00a8e1'
    }

    edge_to_rich_sentence = {}
    # ThemeOf relation between gene and mutations from dbSNP
    for norm_edge, (source, target, evidence) in edge_to_selected_evidence.items():

        norm_source, relation, norm_target = norm_edge

        pmid, sent_id, sentence, tagger_word = evidence

        source_type = id_to_type(norm_source)
        target_type = id_to_type(norm_target)

        source_mention, _, source_offset, _ = source
        target_mention, _, target_offset, _ = target

        # (symbol, 'Gene', 'dbSNP', f'{symbol}-Gene-{cor_entrez}')
        if source_offset == 'dbSNP':
            edge_to_rich_sentence[ norm_edge ] = 'Gene-mutation associations from dbSNPs.'
            continue

        # source_start, source_end = eval(source_offset)
        # target_start, target_end = eval(target_offset)

        source_offset = eval(source_offset)
        target_offset = eval(target_offset)

        source_color = type_to_color[source_type]
        target_color = type_to_color[target_type]

        # fixme:
        if relation == 'ThemeOf':
            offset_to_color = {
                source_offset: source_color,
                target_offset: target_color
            }
        else:
            tagger, tagger_offset = tagger_word
            # tagger_start, tagger_end = eval(tagger_offset)
            tagger_offset = eval(tagger_offset)
            tagger_color = type_to_color[relation]

            offset_to_color = {
                source_offset: source_color,
                target_offset: target_color,
                tagger_offset: tagger_color
            }

        sorted_offset = sorted(offset_to_color.keys(), key=lambda x: x[0], reverse=True)
        character_list = [s for s in sentence]

        # start from last offset
        for offset in sorted_offset:
            start, end = offset
            color = offset_to_color[offset]

            character_list.insert(end, '</b></font>')
            character_list.insert(start, f'<font color="white" style="background:{color}"><b>')

        rich_text = ''.join(character_list)

        # 20220622
        # --> conflict with <!-- -->
        rich_text = rich_text.replace('-->', '--\>')

        edge_to_rich_sentence[norm_edge] = (pmid, sent_id, rich_text)

    return edge_to_rich_sentence



def edge_to_echarts(edge_count: dict, norm_id_to_name: dict,
                    edge_to_rich_evidence: dict,
                    only_option: bool=False):
    relation_to_color = {
        # grey
        'ThemeOf': '#666666',
        # black
        'Reg': '#0F0F0F',
        # red
        'PosReg': '#FF6666',
        # blue
        'NegReg': '#00a8e1',
    }

    relation_to_arrow = {
        'ThemeOf': [ '', '' ],
        'Reg': [ '', 'triangle' ],
        'PosReg': [ '', 'triangle' ],
        'NegReg': [ '', 'triangle' ],
    }

    # relation_to_value = {
    #     'ThemeOf': 50,
    #     'Reg': 1,
    #     'PosReg': 1,
    #     'NegReg': 1,
    # }

    edges = []
    for edge in edge_count.keys():

        # todo: edge from dbSNP
        if len(edge_to_rich_evidence[edge]) == 3:
            pmid, sent_id, rich_evidence = edge_to_rich_evidence[edge]
            value = f'{pmid}: {rich_evidence}'
        else:
            value = f'{edge_to_rich_evidence[edge]}'

        source, relation, target = edge
        # print(source)
        if norm_id_to_name.get(source):
            source_name = norm_id_to_name[source]
        else:
            print(f'wrong source: {edge}')
            continue

        if norm_id_to_name.get(target):
            target_name = norm_id_to_name[target]
        else:
            print(f'wrong target: {edge}')
            continue

        color = relation_to_color[ relation ]
        arrow_type = relation_to_arrow[relation]
        # value = relation_to_value[relation]

        # curveness = relation_to_curveness[relation]

        edges.append({
            'source': source_name,
            'target': target_name,
            'value': value,
            # 'value': f'{source_name}-{relation}-{target_name}',
            'symbol': arrow_type,
            'symbolSize': [0, 20],
            'lineStyle': {'color': color, 'curveness': 0.5, 'width': 5},
            'emphasis': {
                'lineStyle': {'width': 10, 'opacity': 0.8},
            },
        })
    if only_option:
        return f"'links': {edges}"
    else:
        return f'var links={edges}'

def node_to_echarts(node_count: dict, norm_id_to_name: dict, only_option:bool=False):

    type_to_node_color = {
        'Gene': '#ffa235',
        'Var': '#8932a5',
        'GO': 'rgba(255,255,255,0)',
        'HPO': 'rgba(255,255,255,0)',
        'Disease': 'rgba(255,255,255,0)'
    }

    type_to_label_backgroundColor = {
        'Gene': '#ffa235',
        'Var': 'rgba(255,255,255,0)',
        'GO': '#90ed7d',
        'HPO': '#70ad47',
        'Disease': '#82B0D2',
    }

    type_to_label_padding = {
        'Gene': 0,
        'Var': 0,
        'GO': [2,3,2,3],
        'HPO': [2,3,2,3],
        'Disease': [2,3,2,3],
    }

    type_to_border_color = {
        'Gene': '#000000',
        'Var': '#000000',
        'GO': '#336633',
        'HPO': '#336633',
        'Disease': '#2878B5',
    }

    type_to_node_border_width = {
        'Gene': 0.5,
        'Var': 0,
        'GO': 0,
        'HPO': 0,
        'Disease': 0,
    }

    type_to_label_border_width = {
        'Gene': 0,
        'Var': 0,
        'GO': 1,
        'HPO': 1,
        'Disease': 1,
    }

    type_to_shape = {
        'Gene': 'roundRect',
        'Var': 'circle',
        'GO': 'circle',
        'HPO': 'circle',
        'Disease': 'circle'
    }

    type_to_label_show = {
        'Gene': 'true',
        'Var':'false',
        'GO': 'true',
        'HPO': 'true',
        'Disease': 'true'
    }

    saved_set = set()

    if only_option:
        data_str = "'data': ["
    else:
        data_str = "var data=["
    for node in node_count.keys():
        if norm_id_to_name.get(node):
            node_name = norm_id_to_name[node]
        else:
            node_name = node

        node_type = id_to_type(node)

        # if 'Gene' in node:
        #     node_type = 'Gene'
        # elif 'Var' in node:
        #     node_type = 'Var'
        # elif 'GO' in node:
        #     node_type = 'GO'
        # elif 'HP' in node:
        #     node_type = 'HPO'
        # else:
        #     node_type = ''

        color = type_to_node_color[node_type]
        shape = type_to_shape[node_type]

        if node_type == 'Var':
            node_size = 10
            # label_size = 10
        elif node_type == 'Gene':
            node_size = [len(node_name)*14, 40]
            # label_size = [len(node_name)*14, 40]
        else:
            node_size = 44
            # label_size = [len(node_name)*14, 40]

        label_show = type_to_label_show[node_type]

        label_background_color = type_to_label_backgroundColor[node_type]
        label_padding = type_to_label_padding[node_type]
        label_border_width = type_to_label_border_width[node_type]

        border_color = type_to_border_color[node_type]

        node_border_width = type_to_node_border_width[node_type]

        if node_name in saved_set:
            continue
        saved_set.add(node_name)

        # 让基因和突变不根据图例取消
        # if node_type == 'Gene':
        #     node_type = 'gene'
        # if node_type == 'Var':
        #     node_type = 'var'

        node_str = "{" + \
            "'name': \"{0}\",\
            'value': '{1}',\
            'symbol': '{2}',\
            'symbolSize': {3},\
            'draggable': 'True',\
            'category': '{4}',\
            'label': ".format(node_name, node, shape, node_size, node_type) + \
            "{" +\
            "'show': {0},\
              'position': 'inside',\
              'fontWeight': 'bolder',\
              'backgroundColor': '{1}',\
              'padding': {2},\
              'borderColor': '{3}',\
              'borderWidth': {4}"\
                       .format(label_show, label_background_color, label_padding, border_color, label_border_width) +\
            "}," +\
            "'itemStyle': {" +\
                "'color': '{0}',\
                 'borderColor': '{1}',\
                 'borderWidth': {2}".format(color, border_color, node_border_width) +\
            "}," \
            "'emphasis': {'label': {'show': true}},\
        }, "

        data_str += node_str

    data_str += "]"

    return data_str


def generate_echarts(node_count: dict, edge_count: dict,
                     edge_to_rich_evidence: dict,
                     temp_path: str='../temp_dir',
                     norm_id_to_name: dict=None,
                     only_option: bool=True,
                     prefix: str=''):

    node_str = node_to_echarts(node_count, norm_id_to_name, only_option)
    edge_str = edge_to_echarts(edge_count, norm_id_to_name, edge_to_rich_evidence, only_option)

    if prefix:
        node_echarts_temp = f'{temp_path}/{prefix}.node_echarts.txt'
        edge_echarts_temp = f'{temp_path}/{prefix}.edge_echarts.txt'
    else:
        node_echarts_temp = f'{temp_path}/node_echarts.txt'
        edge_echarts_temp = f'{temp_path}/edge_echarts.txt'

    with open(node_echarts_temp, 'w') as wf:
        wf.write(f'{node_str}\n')

    with open(edge_echarts_temp, 'w') as wf:
        wf.write(f'{edge_str}\n')

    return node_str, edge_str

def write_html(node_str: str, edge_str: str, save_html: str):
    html_head = """
            <!DOCTYPE html>
            <html lang="en">
            <head>
            <meta charset="UTF-8">
            <script type="text/javascript" src="https://assets.pyecharts.org/assets/echarts.min.js"></script>
            <style>
                .container {
                    margin-left:-100px;
                }
                .container2 {
                    margin-left: 100px;
                    {#margin-right: -300px;#}
                    margin-right: 100px;
                }
                thead {
                    background-color: rgba(244, 67, 54, 0.85);
                    color: #fff;
                }
                mark {
                    background-color: #f4f57a;
                    padding: 0;
                }
                button.dt-button, div.dt-button, a.dt-button {
                    padding: 0.3em 0.3em;
                }
            </style>
            <title>Beautiful Tables</title>
            </head>
            <body>
            
            <center>
        
            </center>
            
            <div id="cc51ebcfed2642fbbc21ce3b95b1f0a7" class="chart-container" style="width: 2400px; height: 2400px;"></div>
            <script>
            function my$(id){
            return document.getElementById(id);
            }
            function figure$(id,data,link){
                var MyChart = echarts.init(document.getElementById(id), 'light', {renderer: 'canvas'});
                var option = {
                        "tooltip": {
                                "show": true,
                                'position': {left: '5%', top: '10%'},
                                "trigger": "item",
                                "triggerOn": "click",
                                "enterable": true,
                                "alwaysShowContent": true,
                                "renderMode": 'html',
                                "width": 20,
                                "overflow": "breakAll",
        
                                formatter: function(param) {
                                    var text = ''
                                    text += param.name
                                    text += '<br/>'
                                    text += param.value
                                    return text
                                },
        
                                "backgroundColor": '#F5F5F5',
                                "borderColor": '#333',
                                "borderWidth": 2,
                                "padding": 10,
        
                                "textStyle":{
                                    'color':"#333",
                                    "fontStyle": 'italic',
                                    "fontWeight": 'normal',
                                    "fontFamily": 'Courier New',
                                    //"fontFamily": 'Arial',
                                    "fontSize": 20,
                                    "width": 30,
                                    "overflow": 'break',
                                    "rich": {
                                    "protein":{"color": "red"},
                                    }
                                }},
                        "textStyle":{
                            "color": "#000000",
                            "fontStyle": 'italic',
                            "fontWeight": 'normal',
                            "fontFamily": 'Courier New',
                            "fontSize": 20,},
        
                        "toolbox": {
                            "show": true,
                            "itemSize": 15,
                            "itemGap": 16,
                            "right": '5%',
                            "top": '5%',
        
                            "feature": {
                                "saveAsImage": {
                                    "show": true,
                                    "type": 'png',
                                    "title": "Save as image.",
                                },
                                "restore": {
                                    "show": true,
                                    "title": "Revert",
                                    },
                                "dataView": {
                                    "show": true,
                                    "title": "Show Data",
                                    "readOnly": true,
                                }}},
        
                "animation": true,
                "animationThreshold": 2000,
                "animationDuration": 1000,
                "animationEasing": "quinticlnOut",
                "animationDelay": 0,
                "animationDurationUpdate": 300,
                "animationEasingUpdate": "cubicOut",
                "animationDelayUpdate": 0,
        
                "color": [
                            "#FCCF73",
                            "#B9D92C",
                            ],
                "series": [
                    {
                        "type": "graph",
                        "layout": "force",
                        "symbolSize": 10,
                        "nodeScaleRatio": 0.6,
                        "circular": {
                            "rotateLabel": true
                        },
                        "force": {
                            "repulsion": 5000,
                            "edgeLength": [1, 50],
                            "gravity": 0.0001,
                            "friction": 0.1
                        },
                        "label": {
                            "show": true,
                            "position": {left: 10, top: '10%'},
                            "margin": 8
                        },
                        "lineStyle": {
                            "show": true,
                            "width": 2,
                            "opacity": 1,
                            "curveness": 0.15,
                            "type": "solid"
                        },
                        "roam": true,
                        "draggable": true,
                        // show adjacency node
                        "focusNodeAdjacency": true,
                        "data": data,
                        "categories": [
                            {
                                "name": "Gene"
                            },
                            {
                                "name": "Phenotype"
                            },
                            {
                                "name": "Protein"
                            },
                            {
                                "name": "Enzyme"
                            },
                            {
                                "name": "Pathway"
                            },
                            {
                                "name": "Interaction"
                            },
                            {
                                "name": "CPA"
                            },
                            {
                                "name": "MPA"
                            },
                            {
                                "name": "Var"
                            }
                        ],
                        "edgeLabel": {
                            "show": false,
                            "position": 'top',
                            "margin": 8
                        },
                        "edgeSymbol": [
                            null,
                            null
                        ],
                        "edgeSymbolSize": [10, 10],
                        "links": link
                    }
                ],
                "legend":
                    {
                        "data": [
                            "Gene",
                            "Phenotype",
        
                        ],
                    "selected": {
                            "CPA": true,
                            "Enzyme": true,
                            "Interaction": true,
                            "MPA": true,
                            "pathway": true,
                            "PosReg": true,
                            "Protein": true,
                            "Gene": true
                        },
                        "show": false,
                        "padding": 5,
                        "itemGap": 20,
                        "itemWidth": 30,
                        "itemHeight": 25,
                        "left": '5%',
                        "top": '5%',
                        "itemStyle": {
                            "borderColor": '#000000',
                            "borderWidth": 1.6,
                        },
                        "textStyle":{
                        'color':"#000000",
                        "fontStyle": 'italic',
                        "fontWeight": 'normal',
                        "fontFamily": 'Courier New',
                        "fontSize": 16,},
                    },
            "title": [
                    {
                        "padding": 5,
                        "itemGap": 10
                    }
                ],
            };
                    MyChart.setOption(option);
        }
    
    """

    html_tail = """
    figure$('cc51ebcfed2642fbbc21ce3b95b1f0a7',data,links);
    </script>
    """

    with open(save_html, 'w')  as wf:
        wf.write(f'{html_head}\n')
        wf.write(f'{node_str}\n')
        wf.write(f'{edge_str}\n')
        wf.write(f'{html_tail}\n')

    print(f'{save_html} save done.')

def read_html(html_frag_file: str):
    html_doc = open(html_frag_file).read()
    return html_doc

def get_table_data(edge_count, norm_id_to_name, edge_to_rich_evidence):

    table_to_count = defaultdict(int)
    for edge in edge_count.keys():
        source, relation, target = edge
        if relation == 'ThemeOf':
            continue
        print(edge)
        norm_source = norm_id_to_name[source]
        entrez = source.split('-')[0]
        symbol, rs_id = norm_source.split(':')

        norm_target = norm_id_to_name[target]
        if '-' in target:
            term_id = target.split('-')[1]
        else:
            term_id = target

        table_vec = (entrez, symbol, rs_id, relation, term_id, norm_target)

        print(table_vec)
        input()
        table_to_count[table_vec] += 1

    return table_to_count

def generate_graph_option(chain_file: str,
                          graph_html_file: str='../data/html_fragments/graph_option.txt',
                          only_option: bool=True):
    """
    data_str_put_here
    edge_str_put_here
    """
    node_count, edge_count, norm_id_to_name, edge_to_rich_evidence = read_chain_file(chain_file)

    table_data = get_table_data(edge_count, norm_id_to_name, edge_to_rich_evidence)


    node_str, edge_str = generate_echarts(node_count=node_count,
                                          edge_count=edge_count,
                                          edge_to_rich_evidence=edge_to_rich_evidence,
                                          norm_id_to_name=norm_id_to_name,
                                          only_option=only_option,
                                          temp_path='../temp_dir',
                                          )

    graph_html = read_html(graph_html_file)

    graph_html = graph_html.replace('data_str_put_here', node_str)

    graph_html = graph_html.replace('edge_str_put_here', edge_str)

    return graph_html



def main():
    parser = argparse.ArgumentParser(description='SemanticRoleChain to Echarts.')

    parser.add_argument('-cf', dest='chain_file', required=True,
                        help='Semantic Role Chain file.')
    parser.add_argument('-sh', dest='save_html', required=True,
                        help='saved html file.')
    parser.add_argument('-td', dest='temp_dir', default='../temp',
                        help='Temporary folder, default: ../temp.')
    parser.add_argument('-pf', dest='prefix', default='',
                        help='save prefix.')


    # parser.add_argument('-gf', dest='go_file',
    #                     default='../data/obo-data/go.obo',
    #                     help='default: ../data/obo-data/go.obo, download it.')
    #
    # parser.add_argument('-if', dest='gene_info_file',
    #                     default='../data/ncbi-data/Homo_sapiens.gene_info.gz',
    #                     help='default: ../data/ncbi-data/Homo_sapiens.gene_info.gz, download it.')

    args = parser.parse_args()

    node_count, edge_count,\
    norm_id_to_name, edge_to_rich_evidence = read_chain_file(args.chain_file)

    node_str, edge_str = generate_echarts(node_count=node_count,
                                          edge_count=edge_count,
                                          edge_to_rich_evidence=edge_to_rich_evidence,
                                          norm_id_to_name=norm_id_to_name,
                                          temp_path=args.temp_dir,
                                          prefix=args.prefix,
                                          )

    # write_html(node_str, edge_str, args.save_html)

if __name__ == '__main__':
    main()
