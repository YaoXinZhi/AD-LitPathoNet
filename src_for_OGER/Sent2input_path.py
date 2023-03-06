# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 16/06/2022 11:46
@Author: XINZHI YAO
"""

import os
import argparse

"""
改代码用于将一个文件分割到多个文件夹
每个文件一句话
一个文件夹特定数量文件
方便愚蠢的手动多进程
"""


def sent_to_input(sent_file: str, save_dir: str, dir_size: int,
                  old_path: str):

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    if old_path:
        exist_file_set = set(os.listdir(old_path))
    else:
        exist_file_set = set()

    dir_count = 0
    dir_num = 0
    save_path = f'{save_dir}/dir_{dir_num}'

    if not os.path.exists(save_path):
        os.mkdir(save_path)
    with open(sent_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')

            sent_id = l[0]
            pmid = l[1]
            sentence = l[2]

            # '35370546.2520461.txt'
            # pmid.sent_id.txt

            file_name = f'{pmid}.{sent_id}.txt'

            # fixme: 20220618
            # generate input file for uncompleted tag file.
            if file_name not in exist_file_set:
            # if file_name in exist_file_set:
                continue

            save_file = f'{save_path}/{file_name}'

            with open(save_file, 'w') as wf:
                wf.write(f'[{sentence}]')

            dir_count += 1

            if dir_count % dir_size == 0:
                print(f'{save_path} save done.')
                dir_num += 1
                save_path = f'{save_dir}/dir_{dir_num}'

                if not os.path.exists(save_path):
                    os.mkdir(save_path)


def main():

    parser = argparse.ArgumentParser('sent 2 input path.')
    parser.add_argument('-if', dest='sent_file', required=True)
    parser.add_argument('-sd', dest='save_dir', required=True)
    parser.add_argument('-ds', dest='dir_size', default=1000000, type=int )
    parser.add_argument('-op', dest='old_path', default='')
    args = parser.parse_args()

    sent_to_input(args.sent_file, args.save_dir, args.dir_size, args.old_path)



if __name__ == '__main__':
    main()
