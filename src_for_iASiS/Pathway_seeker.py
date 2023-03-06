# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 22/09/2022 21:25
@Author: XINZHI YAO
"""

import re
import os
import argparse
import networkx as nx
from collections import defaultdict

class concept:
    def __init__(self, node_label: str, node_type:str, node_id: str):

        self.node_label = node_label
        self.node_type = node_type
        self.node_id = node_id

def nx_test():

    # Networkx document: https://networkx.org/documentation/stable/index.html
    DG = nx.DiGraph()

    # add node with attribution
    DG.add_nodes_from([(3, {'type': 'Gene'})])

    # G.add_edges_from([(2, 4, {'weight': 3})])
    DG.add_edges_from([(1,2, {'type': 'PosReg'}), (3,1, {'type': 'EFFECT'})])

    # check the node information
    DG.nodes().data()
    # check the edge information
    DG.edges()
    # all adj-node and edge attribution of node 1
    print(DG[1])

    print(DG[1][2])

    # all edge with attribution
    for (u, v, att) in DG.edges.data():
        print(u, v, att)

    # draw the graph
    nx.draw_networkx(DG)
    # plt.show()

    # shortest path
    sp = nx.shortest_path(DG, 3, 2)
    # Create a graph from 'sp'
    pathGraph = nx.path_graph(sp)  # does not pass edges attributes

    for ea in pathGraph.edges():
        # print from_node, to_node, edge's attributes
        print(ea, DG.edges[ ea[ 0 ], ea[ 1 ] ])

class PNRLE_network:
    def __init__(self, chain_file: str, iasis_file: str,
                 save_path: str, highest_level: str,
                 save_count: int, add_iasis_mode: str,
                 prefix: str, no_norm_var_file: str,
                 rs_entrez_symbol_file: str,
                 add_non_norm_var: bool):

        self.chain_file = chain_file
        self.iasis_file = iasis_file
        self.save_path = save_path
        self.prefix = prefix

        self.highest_level = highest_level
        self.save_count = save_count

        self.add_iasis_mode = add_iasis_mode

        self.add_non_norm_var = add_non_norm_var
        self.no_norm_var_file = no_norm_var_file
        self.rs_entrez_symbol_file = rs_entrez_symbol_file

        # concept_id 2 class node
        self.id_to_node = {}
        # triple set used to initialize the graph
        self.triple_set = set()
        # concept_id set use to filter the iasis data
        self.concept_set = set()

        # statistic of relation type
        self.relation_count = defaultdict(int)

        self.gene_set = set()
        self.go_set = set()
        self.hpo_set = set()
        self.disease_set = set()

        self.no_norm_var_set = set()
        self.rs_to_entrez_symbol = {}

        self.type_order = [ 'GO', 'HPO', 'Disease' ]

        self.ban_type = []
        self.save_type = []

        self.DG = nx.DiGraph()

        self.type_filter__init__()

        if not os.path.exists(save_path):
            os.mkdir(save_path)

    def read_rs_entrez_symbol_file(self,):

        with open(self.rs_entrez_symbol_file) as f:
            f.readline()
            for line in f:
                rs, entrez, symbol = line.strip().split('\t')
                self.rs_to_entrez_symbol[ rs ] = (entrez, symbol)

    def read_non_var_file(self):
        with open(self.no_norm_var_file) as f:
            f.readline()
            for line in f:
                var_list = line.strip().split('\t')[-1].split('|')
                self.no_norm_var_set.update(var_list)

    def type_filter__init__(self):
        self.ban_type = self.type_order[self.type_order.index(self.highest_level) + 1: ]
        self.save_type = self.type_order[:self.type_order.index(self.highest_level) + 1 ]

    def concept_to_id(self, concept_norm: str, concept_type: str, concept_label: str):

        # print(f'concept norm: {concept_norm}')
        norm_split = concept_norm.split('-')
        concept_id = norm_split[-1]
        concept_name = '-'.join(norm_split[:-2])

        # concept_name, _, concept_id = concept_norm.split('-')
        # print(concept_norm, concept_type)
        if concept_type == 'Disease' and concept_id.startswith('MESH'):
            concept_id = concept_id.split(':')[1]
            concept_label = concept_name

        if concept_type== 'Gene':
            self.gene_set.add(concept_id)
            concept_label = concept_name

        if concept_id.startswith('GO:'):
            self.go_set.add(concept_id)
            concept_label = concept_name

        elif concept_id.startswith('HP:'):
            self.hpo_set.add(concept_id)
            concept_label = concept_name

        elif concept_id.startswith('C') or concept_id.startswith('D'):
            self.disease_set.add(concept_id)
            concept_label = concept_name

        return concept_id, concept_label

    def id_to_type(self, concept_id: str, concept_type: str):
        if concept_type == 'Gene':
            self.gene_set.add(concept_id)
            return 'Gene'
        elif concept_type == 'Var':
            return 'Var'
        elif concept_type in {'Disease', 'Mesh'} \
                and (concept_id.startswith('D') or concept_id.startswith('C')):
            self.disease_set.add(concept_id)
            return 'Disease'
        if concept_id.startswith('GO:'):
            self.go_set.add(concept_id)
            return 'GO'
        elif concept_id.startswith('HP'):
            self.hpo_set.add(concept_id)
            return 'HPO'

        return concept_type

    def triple_update(self, triple: tuple, theme_gene:str=''):
        # print('triple_update')
        # print(triple)
        # input()

        source, relation, target = triple

        source_label, source_type, _, source_norm = source
        target_label, target_type, _, target_norm = target

        if source_type == 'Gene' and re.findall(f'rs\d+', target_norm)\
            and target_type == 'Var' and relation == 'ThemeOf':

            rs_id = re.findall(f'rs\d+', target_norm)[0]
            if not self.rs_to_entrez_symbol.get(rs_id):
                print(f'Missed RsID: {rs_id}')
                return ''

            entrez, symbol = self.rs_to_entrez_symbol[rs_id]

            # _symbol, _, _entrez = source_norm.split('-')
            norm_split = source_norm.split('-')
            _entrez = norm_split[-1]
            _symbol = '-'.join(norm_split[:-2])
            if entrez != _entrez:
                source = (symbol, 'Gene', 'dbSNP', f'{symbol}-Gene-{entrez}')
                source_label, source_type, _, source_norm = source

        if source_type in {'Gene', 'Protein'} and source_norm == '-' \
                and re.findall(f'rs\d+', target_norm) \
                and target_type == 'Var' \
                and relation == 'ThemeOf':
            rs_id = re.findall(f'rs\d+', target_norm)[0]
            if self.rs_to_entrez_symbol.get(rs_id):
                entrez, symbol = self.rs_to_entrez_symbol[rs_id]
            else:
                return False
            source = (symbol, 'Gene', 'dbSNP', f'{symbol}-Gene-{entrez}')
            source_label, source_type, _, source_norm = source

        if source_norm == '-' or target_norm == '-':
            # print('jump out')
            # input()
            return ''

        # fixme: normalized concept name replace concept label
        source_id, source_label = self.concept_to_id(source_norm, source_type, source_label)
        target_id, target_label = self.concept_to_id(target_norm, target_type, target_label)

        # source_id, _ = self.concept_to_id(source_norm, source_type, source_label)
        # target_id, _ = self.concept_to_id(target_norm, target_type, target_label)
        # print(source_id, source_type)
        # print(target_id, target_type)
        # input()

        source_type = self.id_to_type(source_id, source_type)
        target_type = self.id_to_type(target_id, target_type)

        if relation == 'ThemeOf' and source_type == 'Gene' and target_type == 'Var':
            target_id = f'{source_id}:Var-{target_id}'

        if relation in {'Reg', 'PosReg', 'NegReg'} and source_type == 'Var' and theme_gene:
            source_id = f'{theme_gene}:Var-{source_id}'


        # print('triple update 3')
        # print(source)
        # print(source_id)
        # print(source_type)
        # print(source_label)
        # print(source_norm)
        #
        # print(target)
        # print(target_id)
        # print(target_type)
        # print(target_label)
        # print(target_norm)
        # input()

        self.id_to_node[source_id] = concept(source_label, source_type, source_id)
        self.id_to_node[target_id] = concept(target_label, target_type, target_id)

        self.triple_set.add((source_id, relation, target_id))

        self.concept_set.add(source_id)
        self.concept_set.add(target_id)

    def triple_extractor(self, chain: list):
        i=0
        reg_source = ''
        reg_type = ''

        # print('chain in triple_extractor')
        # print(chain)
        # print()
        source_gene_id = ''
        while i+2 < len(chain):
            source = eval(chain[i])
            relation = chain[i+1]
            target = eval(chain[i+2])
            # print('triple extractor loop')
            # print(source)
            # print(relation)
            # print(target)
            # input()

            i += 2
            if source[1] == 'Gene':
                source_gene_id = source[-1].split('-')[-1]

            if relation == 'CauseOf':
                reg_source = source
                reg_type = target[1]
                continue
            # left of reg
            if not reg_source:
                triple = (source, relation, target)
                self.triple_update(triple)
            # right of reg
            else:
                triple = (reg_source, reg_type, target)
                self.triple_update(triple, source_gene_id)
                reg_source = ''
        # print(f'Triple set: {len(self.triple_set):,}')
        # print(f'Concept set: {len(self.concept_set):,}')

    def chain_filter(self, chain: list):

        if not chain:
            return False

        # do not use the variation without normalization.
        if not self.add_non_norm_var:
            # print('no non_var_added')
            if not re.findall(r'rs\d+', str(chain)):
                return False

        # print(chain)

        # gene have to normalized
        gene = eval(chain[0])
        if gene[-1] == '-':
            return False

        # concept have to normalized
        _concept = eval(chain[-1])
        # print('concept no norm filter.')
        if _concept[-1] == '-':
            return False

        # have theme of gene for last concept
        if len(chain) > 5:

            # print('theme filter')
            if chain[-2] == 'ThemeOf' and eval(chain[-1])[1]=='Gene':
                return False

            # ThemeOf 在链条抽取中需要特殊处理 (暂时放弃对复杂PNRLE链条的处理，包括多个ThemeOf等)
            # gene-var 需要看作整体
            # 有些 go themeOf HPO/Disease的不能作为通路进一步连接的依据
            # print('cause filter')
            if chain[-4] != 'CauseOf':
                return False

        if eval(chain[2])[1]=='Var':

            var_norm = eval(chain[2])[-1]

            # only rs id var_norm_id like this
            var_norm_sub = var_norm.split('-')[-1]
            # var_name = eval(chain[2])[0]

            # not rs mutation and not in no_norm_var
            #  wrong normalization for non-norm-var
            # ('polymorphisms', 'Var', '(59, 72)', 'point mutations-Mutation-Mutations')
            # have fixed
            if not re.findall(r'rs\d+', var_norm) \
                    and var_norm_sub.lower() not in self.no_norm_var_set:
                    # and var_name not in self.no_norm_var_set:
                return False
        else:
            # print(f'2 is not var : {eval(chain[2])[1]}')
            return False

        # concept level filter
        for element in chain:
            if element in {'ThemeOf', 'CauseOf'}:
                continue
            if eval(element)[1] in self.ban_type:
                if eval(element)[1] == 'Disease':
                    if re.findall(r'MESH:.*?', element):
                        return False
                    else:
                        return True
                else:
                    return False

        # remove ad in pnrle
        if 'D000544' in str(chain) \
            or 'HP:0002511' in str(chain):
            return False
        return True

    def print_downstream_node(self, source: str):

        for node_id, attr in self.DG[source].items():
            relation = attr['relation']
            node_attr = self.DG.nodes()[node_id]
            node_label = node_attr['label']
            print(f'{relation}--{node_id}--{node_label}')
            input()

    def chain_parser(self):
        with open(self.chain_file) as f:
            for line in f:
                l = line.strip().split('\t')
                # sent_id = l[0]
                chain = l[3:]
                # try:
                if not self.chain_filter(chain):
                    # print(chain)
                    # print('drop out.')
                    # input()
                    continue
                # print(chain)
                # print('keeped')
                # input()
                self.triple_extractor(chain)
                # except Exception as e:
                #     print(f'Error detail:\n {e.__class__.__name__}, {e}')
                #     input()
                #     # print('错误明细是', , e)
                #     continue


    def pnrle_network_init(self):

        for concept_id, node in self.id_to_node.items():

            node_label = node.node_label
            node_type = node.node_type
            node_id = node.node_id
            if node_type == 'Gene':
                if len(node_label) > 10:
                    print(node_id)
                    print(node_label)
                    input()
            # create nodes
            self.DG.add_nodes_from([(concept_id,  {'label': node_label,
                                              'type': node_type,
                                              'id': node_id})])

        for source_id, relation, target_id in self.triple_set:

            self.DG.add_edges_from([(source_id, target_id, {'relation': relation})])
            self.relation_count[relation] += 1

        # print(f'Network initialized, including {len(self.DG.nodes):,} nodes, '
        #       f'{len(self.DG.edges):,} edges.')
        self.print_network_statistic()

    def print_network_statistic(self):
        print(f'Network includes {len(self.DG.nodes):,} nodes; '
              f'{len(self.DG.edges):,} edges.')
        print(f'Gene: {len(self.gene_set):,};')
        print(f'GO: {len(self.go_set):,}; '
              f'HPO: {len(self.hpo_set):,};'
              f'Disease: {len(self.disease_set):,}.')

        relation_statistic = ''
        for relation, count in self.relation_count.items():
            relation_statistic += f'{relation}: {count:,}; '

        print(relation_statistic)

    def iasis_parser(self):

        iasis_add_count = 0
        with open(self.iasis_file) as f:
            f.readline()
            for line in f:
                # print(line.strip().split('\t'))
                # print(len(line.strip().split('\t')))
                if len(line.strip().split('\t')) == 10:
                    relation, relation_element, source_label, source_id, source_type, _, target_label, target_id, target_type, _ = line.strip().split('\t')
                else:
                    relation, source_label, source_id, source_type, _, target_label, target_id, target_type, _ = line.strip().split('\t')

                # do not add new node from iasis
                if self.add_iasis_mode == 'two':
                    if source_id not in self.concept_set \
                            or target_id not in self.concept_set:
                        continue
                # at least one node in pnrle
                elif self.add_iasis_mode == 'one':
                    if source_id not in self.concept_set \
                        and target_id not in self.concept_set:
                        continue
                else:
                    pass
                # at least one node in pnrle
                # if source_id not in self.concept_set \
                #         and target_id not in self.concept_set:
                #     continue

                    # merge alzheimer hpo-mesh
                if source_id == 'HP:0002511':
                    source_id = 'D000544'
                    source_type = 'Disease'
                if target_id == 'HP:0002511':
                    target_id = 'D000544'
                    target_type = 'D000544'

                source_type = self.id_to_type(source_id, source_type)
                self.DG.add_nodes_from([(source_id, {'label': source_label,
                                                    'type': source_type,
                                                    'id': source_id})])
                self.DG.add_nodes_from([(target_id, {'label': target_label,
                                                     'type': target_type,
                                                     'id': target_id})])

                # add edge from iasis
                self.DG.add_edges_from([(source_id, target_id, {'relation': relation})])
                self.relation_count[relation] += 1

                iasis_add_count += 1

        print(f'{iasis_add_count:,} edges from iasis are added.')
        # print(f'Network including {len(self.DG.nodes):,} nodes, '
        #       f'{len(self.DG.edges):,} edges.')
        self.print_network_statistic()

    def print_simple_path(self, source: str, target: str, cutoff: int=float('inf')):
        for path in nx.all_simple_paths(self.DG, source, target, cutoff=cutoff):
            label_path = [ ]
            type_path = [ ]
            id_path = [ ]

            for idx, node in enumerate(path):
                label_path.append(self.DG.nodes[ node ][ 'label' ])
                type_path.append(self.DG.nodes[ node ][ 'type' ])
                id_path.append(self.DG.nodes[ node ][ 'id' ])

                # relation
                if idx + 1 < len(path):
                    relation = self.DG[ node ][ path[ idx + 1 ] ][ 'relation' ]
                    label_path.append(relation)
                    type_path.append('-')
                    id_path.append('-')

            label_wf = '\t'.join(label_path)
            type_wf = '\t'.join(type_path)
            id_wf = '\t'.join(id_path)
            print(label_wf)
            print(type_wf)
            print(id_wf)
            print()
            input('typing for continue.')

    def save_fixed_pathway(self, source: str, target: str,
                           save_file: str='', cutoff: int=float('inf')):
        if not nx.has_path(self.DG, source, target):
            print(f'No path between {source} and {target}.')

        path_generator = nx.all_simple_paths(self.DG, source, target, cutoff=cutoff)

        save_count = 0
        with open(save_file, 'w') as wf:
            for path in path_generator:
                label_path = [ ]
                type_path = [ ]
                id_path = [ ]

                for idx, node in enumerate(path):
                    label_path.append(self.DG.nodes[ node ][ 'label' ])
                    type_path.append(self.DG.nodes[ node ][ 'type' ])
                    id_path.append(self.DG.nodes[ node ][ 'id' ])

                    # relation
                    if idx + 1 < len(path):
                        relation = self.DG[ node ][ path[ idx + 1 ] ][ 'relation' ]
                        label_path.append(relation)
                        type_path.append('-')
                        id_path.append('-')

                path_length = len(path)
                label_wf = '\t'.join(label_path)
                type_wf = '\t'.join(type_path)
                id_wf = '\t'.join(id_path)
                save_count += 1
                wf.write(f'{path_length}\t{source}\t{target}\t{label_wf}\n')
                wf.write(f'{path_length}\t{source}\t{target}\t{type_wf}\n')
                wf.write(f'{path_length}\t{source}\t{target}\t{id_wf}\n')
                wf.write('\n')
                if save_count >= self.save_count:
                    return ''

        print(f'{save_file} save done, {save_count} path saved.')

    def save_pathway(self, source_set: set, target_set: set, save_file: str,
                     cutoff: int):
        save_count = 0
        with open(save_file,  'w') as wf:
            wf.write(f'PathLength\tSource\tTarget\n')
            for source in source_set:
                for target in target_set:
                    if not self.DG.has_node(source) or not self.DG.has_node(target):
                        continue

                    if not nx.has_path(self.DG, source, target):
                        continue

                    for path in nx.all_simple_paths(self.DG, source, target, cutoff=cutoff):

                        label_path = [ ]
                        type_path = [ ]
                        id_path = [ ]

                        for idx, node in enumerate(path):
                            label_path.append(self.DG.nodes[ node ][ 'label' ])
                            type_path.append(self.DG.nodes[ node ][ 'type' ])
                            id_path.append(self.DG.nodes[ node ][ 'id' ])

                            # relation
                            if idx + 1 < len(path):
                                relation = self.DG[node][path[idx + 1]][ 'relation' ]
                                label_path.append(relation)
                                type_path.append('-')
                                id_path.append('-')

                        path_length = len(path)
                        label_wf = '\t'.join(label_path)
                        type_wf = '\t'.join(type_path)
                        id_wf = '\t'.join(id_path)
                        save_count += 1
                        wf.write(f'{path_length}\t{source}\t{target}\t{label_wf}\n')
                        wf.write(f'{path_length}\t{source}\t{target}\t{type_wf}\n')
                        wf.write(f'{path_length}\t{source}\t{target}\t{id_wf}\n')
                        wf.write('\n')
                        if save_count >= self.save_count:
                            return ''
        print(f'{save_file} save done, {save_count:,} pathway saved.')

    def save_shortest_path(self, source_set: set, target_set: set,
                           save_file: str):
        """
        for i in nx.shortest_simple_paths(DG, '5663', 'D000544'):
            print(i)
        """
        save_count = 0
        with open(save_file, 'w') as wf:
            wf.write('PathLength\tSource\tTarget\n')
            for source in source_set:
                for target in target_set:
                    if not self.DG.has_node(source) or not self.DG.has_node(target):
                        continue

                    if not nx.has_path(self.DG, source, target):
                        continue
                    for path in nx.shortest_simple_paths(self.DG, source, target):
                        label_path = [ ]
                        type_path = [ ]
                        id_path = [ ]
                        for idx, node in enumerate(path):
                            label_path.append(self.DG.nodes[ node ][ 'label' ])
                            type_path.append(self.DG.nodes[ node ][ 'type' ])
                            id_path.append(self.DG.nodes[ node ][ 'id' ])
                            # relation
                            if idx + 1 < len(path):
                                relation = self.DG[node][path[idx + 1]][ 'relation' ]
                                label_path.append(relation)
                                type_path.append('-')
                                id_path.append('-')

                        path_length = len(path)
                        label_wf = '\t'.join(label_path)
                        type_wf = '\t'.join(type_path)
                        id_wf = '\t'.join(id_path)
                        save_count += 1
                        wf.write(f'{path_length}\t{source}\t{target}\t{label_wf}\n')
                        wf.write(f'{path_length}\t{source}\t{target}\t{type_wf}\n')
                        wf.write(f'{path_length}\t{source}\t{target}\t{id_wf}\n')
                        wf.write('\n')
                        if save_count >= self.save_count:
                            return ''
        print(f'{save_file} save done, {save_count} pathway saved.')


    def pathway_search(self, cutoff: int, restrict_concept_set: set):
        print('Search simple pathway.')

        if self.prefix:
            gene_go_save_file = f'{self.save_path}/{self.prefix}-Gene_GO.pathway.tsv'
            gene_hpo_save_file = f'{self.save_path}/{self.prefix}-Gene_HPO.pathway.tsv'
            gene_dis_save_file = f'{self.save_path}/{self.prefix}-Gene_Dis.pathway.tsv'
            gene_restrict_concept_file = f'{self.save_path}/{self.prefix}-Gene_RestrictedConcept.pathway.tsv'
        else:
            gene_go_save_file = f'{self.save_path}/Gene_GO.pathway.tsv'
            gene_hpo_save_file = f'{self.save_path}/Gene_HPO.pathway.tsv'
            gene_dis_save_file = f'{self.save_path}/Gene_Dis.pathway.tsv'
            gene_restrict_concept_file = f'{self.save_path}/Gene_RestrictedConcept.pathway.tsv'

        if restrict_concept_set:
            print(f'Saving Gene--RestrictedConcept pathway.')
            self.save_pathway(self.gene_set, restrict_concept_set, gene_restrict_concept_file, cutoff)
            # 2022-10-18 for more pathway recall
            # print(f'Saving Gene--HPO pathway.')
            # self.save_pathway(self.gene_set, self.hpo_set, gene_hpo_save_file, cutoff)
            # print(f'Saving Gene--Dis pathway.')
            # self.save_pathway(self.gene_set, self.disease_set, gene_dis_save_file, cutoff)

            return ''

        if self.highest_level == 'GO':
            print(f'Saving Gene--GO pathway.')
            self.save_pathway(self.gene_set, self.go_set, gene_go_save_file, cutoff)
        elif self.highest_level == 'HPO':
            print(f'Saving Gene--HPO pathway.')
            self.save_pathway(self.gene_set, self.hpo_set, gene_hpo_save_file, cutoff)
        elif self.highest_level == 'Disease':
            print(f'Saving Gene--Dis pathway.')
            self.save_pathway(self.gene_set, self.disease_set, gene_dis_save_file, cutoff)

    def shortest_path_search(self):

        if self.prefix:
            gene_go_shortest_file = f'{self.save_path}/{self.prefix}-Gene_GO.shortest.tsv'
            gene_hpo_shortest_file = f'{self.save_path}/{self.prefix}-Gene_HPO.shortest.tsv'
            gene_Dis_shortest_file = f'{self.save_path}/{self.prefix}-Gene_Dis.shortest.tsv'
        else:
            gene_go_shortest_file = f'{self.save_path}/Gene_GO.shortest.tsv'
            gene_hpo_shortest_file = f'{self.save_path}/Gene_HPO.shortest.tsv'
            gene_Dis_shortest_file = f'{self.save_path}/Gene_Dis.shortest.tsv'

        if self.highest_level == 'GO':
            print('Saving Gene-GO shortest path.')
            self.save_shortest_path(self.gene_set, self.go_set, gene_go_shortest_file)
        elif self.highest_level == 'HPO':
            print('Saving Gene-HPO shortest path.')
            self.save_shortest_path(self.gene_set, self.hpo_set, gene_hpo_shortest_file)
        elif self.highest_level == 'Disease':
            print('Saving Gene-Dis shortest path.')
            self.save_shortest_path(self.gene_set, self.disease_set, gene_Dis_shortest_file)

    def has_link(self, source: str, target: str):

        print(f'{source}->{target}: {(source, target) in self.DG.edges()}')
        print(f'{source}->{target}: {(target, source) in self.DG.edges()}')
        return  (source, target) in self.DG.edges()

    @staticmethod
    def read_restrict_concept_file(restrict_concept_file: str):

        restrict_id_set = set()
        with open(restrict_concept_file) as f:
            for line in f:
                concept_id, name = line.strip().split('\t')
                restrict_id_set.add(concept_id)

        print(f'Restricted concept: {len(restrict_id_set):,}.')
        return restrict_id_set


    def save_root_leaf_link(self):
        if self.prefix:
            save_file = f'{self.save_path}/{self.prefix}-root_leaf.tsv'
        else:
            save_file = f'{self.save_path}/root_leaf.tsv'

        roots = (v for v, d in self.DG.in_degree if d == 0)
        leaves = [v for v, d in self.DG.out_degree if d == 0]

        save_count = 0
        wf = open(save_file, 'w')
        for source in roots:
            for target in leaves:
                if not self.DG.has_node(source) or not self.DG.has_node(target):
                    continue
                if not nx.has_path(self.DG, source, target):
                    continue
                for path in nx.shortest_simple_paths(self.DG, source, target):
                    label_path = [ ]
                    type_path = [ ]
                    id_path = [ ]
                    for idx, node in enumerate(path):
                        label_path.append(self.DG.nodes[ node ][ 'label' ])
                        type_path.append(self.DG.nodes[ node ][ 'type' ])
                        id_path.append(self.DG.nodes[ node ][ 'id' ])
                        # relation
                        if idx + 1 < len(path):
                            relation = self.DG[ node ][ path[ idx + 1 ] ][ 'relation' ]
                            label_path.append(relation)
                            type_path.append('-')
                            id_path.append('-')

                    path_length = len(path)
                    label_wf = '\t'.join(label_path)
                    type_wf = '\t'.join(type_path)
                    id_wf = '\t'.join(id_path)
                    save_count += 1
                    wf.write(f'{path_length}\t{source}\t{target}\t{label_wf}\n')
                    wf.write(f'{path_length}\t{source}\t{target}\t{type_wf}\n')
                    wf.write(f'{path_length}\t{source}\t{target}\t{id_wf}\n')
                    wf.write('\n')
                    if save_count >= self.save_count:
                        return ''
        wf.close()
        print(f'Root-Leave path save done, {save_count} links saved.')


def main(chain_file: str, iasis_file: str, save_path: str,
         cutoff: int, highest_level: str, save_count: int,
         add_iasis_mode: str, restrict_concept_file: str,
         prefix: str, no_norm_var_file: str, rs_entrez_symbol_file: str,
         add_non_norm_var: bool):

    pathway_seeker = PNRLE_network(chain_file, iasis_file,
                                   save_path, highest_level,
                                   save_count, add_iasis_mode,
                                   prefix, no_norm_var_file,
                                   rs_entrez_symbol_file,
                                   add_non_norm_var)
    # read no_norm_var
    pathway_seeker.read_non_var_file()

    pathway_seeker.read_rs_entrez_symbol_file()

    print('Parsing Chain file.')
    pathway_seeker.chain_parser()

    print('Initializing the PNRLE network.')
    pathway_seeker.pnrle_network_init()

    print('Initializing iASiS network.')
    pathway_seeker.iasis_parser()

    if restrict_concept_file:
        restrict_concept_set = pathway_seeker.read_restrict_concept_file(restrict_concept_file)
    else:
        restrict_concept_set = set()

    # cutoff
    # cutoff = 3
    pathway_seeker.pathway_search(cutoff, restrict_concept_set)
    #
    # pathway_seeker.shortest_path_search()

    # pathway_seeker.save_root_leaf_link()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-cf', dest='chain_file', required=True)
    parser.add_argument('-if', dest='iasis_file', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    parser.add_argument('-pr', dest='prefix', default='',
                        help='save prefix, default: "".')


    parser.add_argument('-hl', dest='highest_level',
                        choices=['GO', 'HPO', 'Disease'],
                        default='Disease',
                        help='default: Disease')

    parser.add_argument('-tf', dest='cutoff', type=int, default=3, help='default: 4')
    parser.add_argument('-sc', dest='save_count', type=int, default=float('inf'),
                        help='pathway save count, default: inf.')
    # parser.add_argument('-op', dest='only_pnrle_nodes', action='store_true',
    #                     default=False, help='do not add new node from iASiS.')
    parser.add_argument('-op', dest='add_iasis_mode', choices=['two', 'one', 'zero'],
                        default='zero', help='At least two/one/zero node from PNRLE.')

    # restricted end node of pathway
    parser.add_argument('-rc', dest='restrict_concept_file', default='../result/End_of_pathway.txt',
                        help='Restricted concept file for pathway generation, default: ../result/End_of_pathway.txt.')

    # non-normalized mutation file
    parser.add_argument('-an', dest='add_non_norm_var', action='store_true',
                        default=False,
                        help='add_non_norm_var, default: False')
    parser.add_argument('-nn', dest='no_norm_var_file',
                        default='../result/MutationTypeDictionary20220711.txt',
                        help='default: ../result/MutationTypeDictionary20220711.txt')

    parser.add_argument('-rg', dest='rs_entrez_symbol_file',
                        default='../result/rs_to_entrez/rs-entrez-symbol.merge.tsv',
                        help='default: ../result/rs_to_entrez/rs-entrez-symbol.merge.tsv')

    args = parser.parse_args()

    """
    # small
    chain_file = 'ad.total.only_theme.chain.6-26.tsv'
    # iasis_file = 'PNRLE_from_iASiS.small.tsv'
    iasis_file = 'PNRLE_from_iASiS.small.FunRel.tsv'
    
    #big
    # chain_file = '/home/xzyao/AD-PNRLE/data/ad-731/ad.chain.OnlyTheme.tsv'
    # iasis_file =  '/home/xzyao/AD-PNRLE/data/ad-731/iASiS_convert/PNRLE_from_iASiS.big.tsv'
    # iasis_file =  '/home/xzyao/AD-PNRLE/data/ad-731/iASiS_convert/PNRLE_from_iASiS.big.FunRel.tsv'
    
    restrict_concept_file = '/hone/xzyao/AD-PNRLE/result/End_of_pathway.txt'
    
    save_path = 'pathway_dir'

    # highest_level = 'GO'
    # highest_level = 'HPO'
    highest_level = 'Disease'
    add_iasis_mode = one

    cutoff = 5
    # save_count = 1000000
    """


    main(args.chain_file, args.iasis_file, args.save_path,
         args.cutoff, args.highest_level, args.save_count,
         args.add_iasis_mode, args.restrict_concept_file,
         args.prefix, args.no_norm_var_file,
         args.rs_entrez_symbol_file, args.add_non_norm_var)

