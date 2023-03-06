# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 16/06/2022 17:30
@Author: XINZHI YAO
"""

import os
import re

from tqdm import tqdm

import argparse

from collections import defaultdict

def read_oger_file_new(oger_file, save_origin):
    tag_set = set()
    go_set = set()
    with open(oger_file) as f:
        for line in f:
            l = line.strip().split('\t')

            tag_type = l[ 1 ]

            start = l[ 2 ]
            end = l[ 3 ]

            mention = l[ 4 ]

            preferred = l[ 5 ]
            entity_id = l[ 6 ]
            origin = l[ 9 ]
            if origin in save_origin:
                go_set.add(entity_id)

                tag_set.add((mention, preferred, tag_type, entity_id, (start, end)))
    return tag_set, go_set

def read_sent_file(sent_file: str):
    print('reading sent file.')
    sent_set = set()
    with open(sent_file) as f:
        for line in f:
            l = line.strip().split('\t')
            sent_id, pmid, sentence = l[ 0 ], l[ 1 ], l[ 2 ]

            sent_set.add((sent_id, pmid, sentence))
    return sent_set


def batch_process(sent_file: str, src_dir: str, save_file: str):
    miss_count = 0

    sent_set = read_sent_file(sent_file)
    save_id_type = False
    print('start processing.')
    processed_count = 0
    with open(save_file, 'w') as wf:
        for (sent_id, pmid, sentence) in tqdm(sent_set):
            tag_file = f'{src_dir}/{pmid}.{sent_id}.txt'

            if os.path.exists(tag_file):
                processed_count += 1

                # local get oger use read_oger_file_new
                tag_result, go_set = read_oger_file_new(tag_file, ['Gene Ontology'])

                try:
                    if not save_id_type:
                        sent_id = sent_id.split(':')[ 1 ]
                except:
                    sent_id = sent_id

                if tag_result:
                    tag = str(tag_result)
                    wf.write(f'{sent_id}\t{pmid}\t{sentence}\t{tag}\n')
            else:
                miss_count += 1

    print(f'{save_file} save done.')


def main():
    parser = argparse.ArgumentParser(description='batch oger result process.')
    parser.add_argument('-ef', dest='sent_file', required=True)
    parser.add_argument('-sd', dest='src_dir', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    batch_process(args.sent_file, args.src_dir, args.save_file)

if __name__ == '__main__':
    main()
