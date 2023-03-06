# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 16/06/2022 9:52
@Author: XINZHI YAO
"""

import os

import argparse

from multiprocessing import Process

import time
from oger.ctrl.router import Router, PipelineServer
# from .OGER_py_Tagger import oger_py_tagger

def oger_py_tagger(src_file: str, save_file: str, pl, process_idx: int):
    doc = pl.load_one(src_file, 'txt')

    pl.process(doc)

    with open(save_file, 'w', encoding='utf-8') as wf:
        pl.write(doc, 'tsv', wf)

# def Multithreaded_call(temp_path: str, function: str,
#                        save_path: str):
def Multithreaded_call(temp_path: str,
                       save_path: str,
                       threading_num: int):

    print(f'Start running.')
    conf = Router(termlist_path='testfiles/go.term.tsv')
    pl = PipelineServer(conf)

    if not os.path.exists(save_path):
        os.mkdir(save_path)

    file_list = os.listdir(temp_path)

    pool_list = []
    threading_count = 0
    for process_idx, _file in enumerate(file_list):
        input_file = f'{temp_path}/{_file}'
        save_file = f'{save_path}/{_file}'

        pool_list.append(f'process_{process_idx}')
        # print(f'process: {process_idx}')
        exec(f'process_{process_idx}=Process(target={oger_py_tagger}, '
             f'args=("{input_file}", "{save_file}", {pl}, {process_idx}))')

        threading_count += 1
        if threading_count % threading_num == 0:
            for pool in pool_list:
                exec(f'{pool}.start()')
                print(f'{pool} running.')
                time.sleep(10)
            pool_list = []


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Multithreaded call.')
    parser.add_argument('-ip', dest='iteration_path', type=str,
                        help='Contains the directory of split files that need to be processed iteratively.')
    parser.add_argument('-sp', dest='save_path', type=str,
                        help='Directory for storing multi-process results.')
    parser.add_argument('-tn', dest='threading_num', default=10)
    # parser.add_argument('-fc', dest='target_function', type=str,
    #                     help='The target function needs to accept three parameters:'
    #                          'split_file: str, save_path: str, process_id: str.')
    args = parser.parse_args()

    # Multithreaded_call(args.iteration_path, args.target_function, args.save_path)
    Multithreaded_call(args.iteration_path, args.save_path, args.threading_num)


if __name__ == '__main__':
    main()
