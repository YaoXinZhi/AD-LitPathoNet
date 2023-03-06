# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 26/09/2022 10:18
@Author: XINZHI YAO
"""

import os
import re
import argparse
from tqdm import tqdm
from collections import defaultdict

"""
改代码用于观测产生的通路
"""

class pathway:
    def __init__(self, path_line: str, type_line: str, label_line: str,
                 score: float):

        self.path_line = path_line
        self.type_line = type_line
        self.label_line = label_line

        self.score = score


class PathwayFilter:
    def __init__(self):
        self.multi_id_mapping_dict = {}
        self.delete_id_set = set()

        # GO:123\tGO:234\tHP:123
        # 用于判断是否保存过这条通路
        # 用于判断是否是冗余通路
        self.saved_path_id_str_set = set()

    def pathway_filter(self, type_line: str, highest_level: str, label_line: str):
        """
        different level in filter
        """
        save_flag = True

        # fixme: do not delete the no_norm var
        # if not re.findall(r'rs\d+', label_line):
        #     return False

        type_list = [_type for _type in type_line.split('\t')[3:] if _type != '-']
        label_list = [label for label in label_line.split('\t')[3:] if label != '-']
        label_set = set(label_list)

        # fixme: delete this rule 1013
        # 1. At least two levels concept.
        # if (len(set(type_list)) == 1 and highest_level != 'GO')\
        #         or (len(type_list) !=3 and highest_level == 'DIS'):
        #     save_flag = False
        #     return save_flag

        type_order = ['GO', 'HPO', 'Disease']

        ban_type = type_order[type_order.index(highest_level)+1:]
        # save_type = type_order[:type_order.index(highest_level)+1]

        # 2. Do not exist higher level concept
        if set(type_list)&set(ban_type):
            save_flag = False
            return save_flag

        # type_idx = defaultdict(int)
        # for idx, _type in enumerate(type_list):
        #     # if _type == 'Disease':
        #     #     _type = 'DIS'
        #     if _type not in save_type:
        #         continue
        #     if idx > type_idx[_type]:
        #         type_idx[_type] = idx

        # fixme: 3. concept level from low to high
        # sort_type = sorted(type_idx.keys(), key=lambda x:type_idx[x])
        # if save_type != sort_type:
        #     save_flag = False

        #4. filter the pathway including selected concepts
        if self.delete_id_set&label_set:
            save_flag = False
            return save_flag

        # cannot include the save concept
        if len(label_list) != len(label_set):
            save_flag = False
            return save_flag

        return save_flag

    @staticmethod
    def id_to_type( _id: str):
        if _id.startswith('GO'):
            return 'GO'
        elif _id.startswith('HP'):
            return 'HPO'
        elif _id.startswith('D') or _id.startswith('C') or _id.lower().startswith('mesh'):
            return 'Dis'
        else:
            return ''

    def score_strict_concept(self, label_line: str, strict_concept_set: set,
                             id_to_label: dict):

        label_split = label_line.split('\t')

        pathway_length = int(label_split[0])
        label_list = [label for label in label_split[2:] if label != '-']

        # print(label_list)
        type_list = [self.id_to_type(label) for label in label_split[2:] if self.id_to_type(label)]
        # print('type_list')
        # print(type_list)

        score = 0

        # encourage more GO
        for _type in type_list:
            if _type == 'GO':
                score += 0.5
            if _type == 'HPO':
                score += 0.25
            if _type == 'Dis':
                score += 0.1

        # do not encourage one more dis
        if type_list.count('Dis') > 1:
            score -= 0.2 * (type_list.count('Dis')-1)

        # first is Mesh
        if type_list[0] == 'Dis':
            score -= 1

        # 1. Include a restricted concept plus 1 point.
        for label in label_list:
            if label in strict_concept_set:
                score += 1

        # 2. One extra point for correct grade order.
        level_list = [id_to_label[label] for label in label_list if id_to_label.get(label)]

        # print(level_list)
        # input()

        level_idx = defaultdict(int)
        for idx, level in enumerate(level_list):
            level_idx[level] = idx

        if len(level_idx) < 2:
            return score

        if 'Molecular level' in level_list and 'Symptom level' in level_list:
            if level_idx['Molecular level'] < level_idx['Symptom level']:
                score += 1
            else:
                score -= 1
        if 'Molecular level' in level_list and 'Disease level' in level_list:
            if level_idx['Molecular level'] < level_idx['Disease level']:
                score += 1
            else:
                score -= 1
        if 'Symptom level' in level_list and 'Disease level' in level_list:
            if level_idx['Symptom level'] < level_idx['Disease level']:
                score += 1
            else:
                score -= 1

        # delete this part, because GO:HP:DIS not always mean level absolutely
        # all is restricted concept
        # more grained concept
        # if len(level_idx) == 3 and score > pre_score:
        #     print(label_list)
        #     print(f'Score before level score: {pre_score}')
        #     print(level_list)
        #     print(level_idx)
        #     print(f'score after level score: {score}')
        #     input()
        score /= (pathway_length-2)
        return score

    @staticmethod
    def save_pathway(id_to_pathway: dict, save_file: str):

        sorted_id = sorted(id_to_pathway.keys(), key=lambda x: id_to_pathway[x].score, reverse=True)

        with open(save_file, 'w') as wf:
            wf.write(f'ChainIdx\tScore\tPathLength\tSource\tTarget\tChain\n')
            for pathway_id in sorted_id:
                path = id_to_pathway[pathway_id]
                path_line = path.path_line
                type_line = path.type_line
                label_line = path.label_line
                score = path.score

                wf.write(f'{pathway_id}\t{score}\t{path_line}\n')
                wf.write(f'{pathway_id}\t{score}\t{type_line}\n')
                wf.write(f'{pathway_id}\t{score}\t{label_line}\n')
                wf.write(f'\n')
        print(f'{save_file} save done, {len(id_to_pathway):,} pathway saved.')

    def save_path_or_not(self, path_line: str, label_line: str):

        # path = path_line.split('\t')
        label_list = label_line.split('\t')
        # path_list = path_line.split('\t')

        # multi convert
        new_label_list = []
        for label in label_list:
            if self.multi_id_mapping_dict.get(label):
                new_label_list.append(self.multi_id_mapping_dict[label])
            else:
                new_label_list.append(label)

        label_str = ','.join([label for label in new_label_list[3:] if label !='-'])
        # saved or not
        if label_str in self.saved_path_id_str_set:
            return False, new_label_list

        # fixme: this part is too slow
        # Redundant pathways
        #for saved_label_str in self.saved_path_id_str_set:

            #if label_str in saved_label_str:
                # print('redundant')
                # print(label_str)
                # print(saved_label_str)
                # input()
                #return False, new_label_list

        self.saved_path_id_str_set.add(label_str)

        # print(self.saved_path_id_str_set)
        # print(path_line)
        # input()

        return True, new_label_list

    def Start_pathway_filter(self, input_file: str, save_file: str, highest_level: str,
                             restrict_concept_file: str):

        if restrict_concept_file:
            print(f'Score pathway using restricted concept.')
            id_to_level, restrict_concept_set = self.read_restrict_concept_file(restrict_concept_file)
        else:
            id_to_level, restrict_concept_set = {}, set()

        print(f'Processing.')
        process_count = 0
        # pathway_id_to_pathway = {}
        pathway_id = 0
        with open(input_file) as f, open(save_file, 'w') as wf:
            wf.write(f'ChainIdx\tScore\tPathLength\tSource\tTarget\tChain\n')
            f.readline()
            while True:
                path_line = f.readline().strip()
                type_line = f.readline().strip()
                label_line = f.readline().strip()

                process_count += 1
                if process_count % 50000 == 0:
                    print(f'{process_count} pathway processed.')

                line = f.readline()
                if not line:
                    break

                if not self.pathway_filter(type_line, highest_level, label_line):
                    continue

                save_path_flag, new_label_list = self.save_path_or_not(path_line, label_line)

                if not save_path_flag:
                    continue

                if restrict_concept_set:
                    score = self.score_strict_concept(label_line, restrict_concept_set,
                                         id_to_level)
                else:
                    score = 'None'

                # print(label_line)
                #
                # print('\t'.join(new_label_list))
                # input()
                # pathway_id_to_pathway[len(pathway_id_to_pathway)] = pathway(path_line, type_line, '\t'.join(new_label_list), score)
                label_line = '\t'.join(new_label_list)
                wf.write(f'{pathway_id}\t{score}\t{path_line}\n')
                wf.write(f'{pathway_id}\t{score}\t{type_line}\n')
                wf.write(f'{pathway_id}\t{score}\t{label_line}\n')
                wf.write(f'\n')
                pathway_id += 1

        print(f'{save_file} save done, {pathway_id} pathway saved.')
        # self.save_pathway(pathway_id_to_pathway, save_file)

    @staticmethod
    def read_restrict_concept_file(restrict_file: str):

        concept_to_level = defaultdict(set)
        restrict_concept_set = set()
        level = ''
        with open(restrict_file) as f:
            for line in f:
                concept_id, concept_name = line.strip().split('\t')
                if concept_id == 'Level':
                    level = concept_name
                    continue
                concept_to_level[concept_id] = level
                restrict_concept_set.add(concept_id)

        return concept_to_level, restrict_concept_set

    def read_multi_id_file(self, multi_id_file: str):
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
                    self.multi_id_mapping_dict[_id] = selected_id

    def read_delete_concept_file(self, delete_concept_file: str):

        with open(delete_concept_file) as f:
            for line in f:
                # print(line.strip().split('\t'))
                _id, concept, count = line.strip().split('\t')
                if self.multi_id_mapping_dict.get(_id):
                    self.delete_id_set.add(self.multi_id_mapping_dict[_id])

                    for or_id, selected_id in self.multi_id_mapping_dict.items():
                        if selected_id == self.multi_id_mapping_dict[_id]:
                            self.delete_id_set.add(or_id)
                else:
                    self.delete_id_set.add(_id)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-if', dest='input_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)

    parser.add_argument('-rf', dest='restrict_concept_file',
                        default='',
                        help='default: None, could be: ../result/ClinicalPhenotype.txt')


    parser.add_argument('-hl', dest='highest_level',
                        choices=['GO', 'HPO', 'Disease'],
                        default='Disease', help='default: Disease')

    parser.add_argument('-dc', dest='delete_concept_file',
                        default='../result/delete_concept.txt',
                        help='default: ../result/delete_concept.txt',
                        required=False)

    parser.add_argument('-mi', dest='multi_id_file', required=False,
                        default='../result/multi_id_concept.tsv',
                        help='default: ../result/multi_id_concept.tsv')

    args = parser.parse_args()

    pathway_filter = PathwayFilter()
    pathway_filter.read_multi_id_file(args.multi_id_file)
    pathway_filter.read_delete_concept_file(args.delete_concept_file)

    pathway_filter.Start_pathway_filter(args.input_file, args.save_file, args.highest_level,
         args.restrict_concept_file, )
