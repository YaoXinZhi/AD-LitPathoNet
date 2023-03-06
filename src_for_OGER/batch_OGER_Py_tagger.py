# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 16/06/2022 14:56
@Author: yao
"""

import os

import argparse

def batch_tagger_commend(src_dir: str, save_file: str,
                         max_threading: int):

    dir_list = os.listdir(src_dir)

    threading_count = 0
    with open(save_file, 'w') as wf:
        for _dir in dir_list:

            wf.write(f'nohup python3 OGER_py_Tagger.py --src_path Sent_split_dir/{_dir} --save_path Sent_split_output_dir/{_dir} &\n')
            threading_count += 1
            if threading_count % max_threading ==0:
                wf.write(f'wait\n')
    print(f'{save_file} save done.')

def main():

    parser = argparse.ArgumentParser('Batch OGER py tagger commend.')
    parser.add_argument('-sd', dest='src_dir', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)
    parser.add_argument('-mt', dest='max_threading', default=3)
    args = parser.parse_args()

    batch_tagger_commend(args.src_dir, args.save_file, args.max_threading)


if __name__ == '__main__':

    main()
