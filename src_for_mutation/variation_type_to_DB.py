# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 24/10/2022 22:29
@Author: XINZHI YAO
"""

def main():

    var_file = f'/home/xzyao/AD-PNRLE/result/rs_to_entrez/rs-entrez-symbol-fxn.tsv'
    db_save_file = f'/home/xzyao/AD-PNRLE/result/rs_to_entrez/rs-type.db.tsv'

    non_coding_var_set = {'intron_variant', '3_prime_UTR_variant',
                          'non_coding_transcript_variant', '5_prime_UTR_variant'}

    with open(var_file) as f, open(db_save_file, 'w') as wf:
        f.readline()
        wf.write('rsID\tFXNClass\tIsNonCodingVar\n')
        for line in f:
            l = line.strip().split('\t')
            rs_id = l[0]
            fax_type = l[3]

            fax_type_set = set(fax_type.split(','))
            #
            # print(fax_type_set)
            # print(fax_type_set & non_coding_var_set)
            # input()

            if fax_type_set&non_coding_var_set:
                have_non_coding_var = 'true'
                # print(rs_id)
                # print(fax_type)
                # print(fax_type_set&non_coding_var_set)
                # print()
            else:
                have_non_coding_var = 'false'

            wf.write(f'{rs_id}\t{fax_type}\t{have_non_coding_var}\n')






if __name__ == '__main__':
    main()
