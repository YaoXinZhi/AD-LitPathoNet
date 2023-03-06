# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 01/07/2022 9:46
@Author: XINZHI YAO
"""

import os
import argparse


def main(sent_file: str, src_path: str, save_path: str, max_threads: int, bash_file: str):
    """
    nohup python3 batch_oger_result_preocess.py
    -ef ../../../../ad-20220415-dir/AD.sent.txt
    -sd Sent_split_output_dir/dir_2
    -sf Sent_split_processed_output/dir_2.tsv &
    """

    src_dir_list = os.listdir(src_path)

    threads_count = 0
    with open(bash_file, 'w') as wf:
        for _dir in src_dir_list:

            src_dir = f'{src_path}/{_dir}'

            save_file = f'{save_path}/{_dir}.tsv'

            if os.path.exists(save_file):
                continue

            commend = f'nohup python3 batch_oger_result_preocess.py -ef {sent_file} -sd {src_dir} -sf {save_file} &'

            wf.write(f'{commend}\n')
            threads_count += 1

            if threads_count % max_threads == 0:
                wf.write(f'wait\n')
    print(f'{bash_file} save done.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='git batch bash result process')
    parser.add_argument('-ef', dest='sent_file',
                        default='../../../../ad-20220415-dir/AD.sent.txt')
    parser.add_argument('-rp', dest='src_path', required=True)
    parser.add_argument('-sp', dest='save_path', required=True)
    parser.add_argument('-mt', dest='max_threads', default=4, type=int)
    parser.add_argument('-sf', dest='bash_file', required=True)
    args = parser.parse_args()

    main(args.sent_file, args.src_path, args.save_path, args.max_threads, args.bash_file)