# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 21/09/2022 17:31
@Author: XINZHI YAO
"""

import os
from collections import defaultdict
import argparse


def main(id_file: str, save_path: str):

	if not os.path.exists(save_path):
		os.mkdir(save_path)

	save_count = 0
	with open(id_file) as f:
		for line in f:
			_id = line.strip()

			save_file = f'{save_path}/{_id}.iASiS.tsv'

			if os.path.exists(save_file):
				# print(f'{save_file} existed.')
				# input()
				continue

			# python iasis_search.py -kw C0000731 -sf ../c0000731.iasis.tsv
			commend = f'python iasis_search.py -kw {_id} -sf {save_file}'

			os.system(commend)
			save_count += 1

			if save_count % 100 == 0:
				print(f'saved concept: {save_count:,}')
			#print(_id)
	print(f'Total {save_count} concepts downloaded.')

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-if', dest='id_file', required=True)
	parser.add_argument('-sp', dest='save_path', required=True)
	args = parser.parse_args()

	"""
	if 
	D:\Mac2PC\AD-PNRLE\data\small_concept_mapped_dir/umls.big.txt
	D:\Mac2PC\AD-PNRLE\data\small_concept_mapped_dir/umls.small.txt
	D:\Mac2PC\AD-PNRLE\data\small_concept_mapped_dir/umls.all.tsv
	sp
	H:\AD_PNRLE_DATA_TEMP\iASiS_all_GMM_download
	"""

	"""
	python Bluk_iASiS_query.py -if D:\Mac2PC\AD-PNRLE\data\small_concept_mapped_dir/umls.all.tsv -sp H:\AD_PNRLE_DATA_TEMP\iASiS_all_GMM_download
	"""


	main(args.id_file, args.save_path)

