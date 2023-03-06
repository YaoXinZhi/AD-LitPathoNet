# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 17/10/2022 20:36
@Author: XINZHI YAO
"""
import argparse
import os
from collections import defaultdict

def read_nc_var_file(nc_var_file: str):

    nc_rs_to_var_type = {}
    with open(nc_var_file) as f:
        f.readline()
        for line in f:
            l = line.strip().split('\t')
            rs_id = l[-1]
            maf_type = l[-3]
            nc_rs_to_var_type[rs_id] = maf_type
    print(f'{len(nc_rs_to_var_type):,} no-coding variation readed.')
    return nc_rs_to_var_type


def main(rs_file: str, nc_var_file: str, save_file: str):

    nc_rs_to_var_type = read_nc_var_file(nc_var_file)

    match_count = 0
    total_count = 0
    with open(rs_file) as f, open(save_file, 'w') as wf:
        head = f.readline().strip()
        wf.write(f'{head}\tMAF_type\n')
        for line in f:
            rs_id, entrez, symbol = line.strip().split('\t')
            total_count += 1
            if nc_rs_to_var_type.get(rs_id):
                match_count += 1
                wf.write(f'{line.strip()}\t{nc_rs_to_var_type[rs_id]}\n')
            else:
                wf.write(f'{line.strip()}\tNone\n')
    print(f'{save_file} saved, {match_count}/{total_count} are non-coding mutations.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-rf', dest='rs_file',
                        default='/home/xzyao/AD-PNRLE/result/rs_to_entrez/rs-entrez-symbol.merge.tsv',
                        help='default: /home/xzyao/AD-PNRLE/result/rs_to_entrez/rs-entrez-symbol.merge.tsv')
    parser.add_argument('-nf', dest='nc_var_file',
                        default='/home/xzyao/AD-PNRLE/result/ncVarDB/ncVarDB-master/data/ncVar_benign.tsv',
                        help='default: /home/xzyao/AD-PNRLE/result/ncVarDB/ncVarDB-master/data/ncVar_benign.tsv')
    parser.add_argument('-sf', dest='save_file',
                        required=True)
    args = parser.parse_args()

    main(args.rs_file, args.nc_var_file, args.save_file)


