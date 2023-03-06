# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 15/04/2022 21:25
@Author: XINZHI YAO
"""

import argparse
import os.path

from oger.ctrl.router import Router, PipelineServer


def oger_py_tagger(src_file: str, save_file: str, pl):
    doc = pl.load_one(src_file, 'txt')

    pl.process(doc)

    with open(save_file, 'w', encoding='utf-8') as wf:
        pl.write(doc, 'tsv', wf)

def batch_tagger(src_path: str, save_path: str, old_path: str=''):
    print(f'Start running.')
    conf = Router(termlist_path='testfiles/go.term.tsv')
    pl = PipelineServer(conf)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    save_count = 0
    for _file in os.listdir(src_path):
        src_file = f'{src_path}/{_file}'
        if old_path:
            old_file = f'{old_path}/{_file}'
        else:
            old_file = ''
        save_file = f'{save_path}/{_file}'
        if os.path.exists(old_file):
            continue
        save_count += 1
        if save_count % 5000 == 0:
            print(save_count)
        oger_py_tagger(src_file, save_file, pl)

    print('Done.')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='OGER python Tagger.')
    parser.add_argument('--src_path', dest='src_path',
                        required=True,
                        help='input file path, each file only include ont sentence with [].')
    parser.add_argument('--save_path', dest='save_path',
                        required=True,
                        help='save path.')
    parser.add_argument('-op', dest='old_path', default='',
                        required=False,
                        help='existed tagging result.')

    args = parser.parse_args()
    batch_tagger(args.src_path, args.save_path, args.old_path)





