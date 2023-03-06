# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 17/06/2022 16:45
@Author: XINZHI YAO
"""

import os

import re
import argparse

from collections import defaultdict

def var_filter(src_file: str, keyword: str,
               count_threshold: int, save_file: str):

    save_count = 0
    with open(src_file) as f, open(save_file, 'w') as wf:
        wf.write(f'Variation\tCount\n')
        f.readline()
        for line in f:
            var, count = line.strip().split('\t')
            count = int(count)

            if count < count_threshold:
                # sorted var file.
                break

            if re.findall(f'{keyword}', var.lower()):
                save_count += 1
                wf.write(f'{var}\t{count}\n')

    print(f'Keyword: {keyword}, save count: {save_count}.')
    print(f'{save_file} save done.')



def main():
    parser = argparse.ArgumentParser(description='No-normalized Variation filter.')
    parser.add_argument('-if', dest='src_file', required=True)
    parser.add_argument('-kw', dest='keyword', required=True)
    parser.add_argument('-ct', dest='count_threshold',
                        default=100, type=int,
                        help='default: 100')
    parser.add_argument('-sf', dest='save_file', required=True)
    args = parser.parse_args()

    var_filter(args.src_file, args.keyword, args.count_threshold, args.save_file)


if __name__ == '__main__':
    main()


