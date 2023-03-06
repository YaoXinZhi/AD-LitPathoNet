# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 24/05/2022 10:54
@Author: XINZHI YAO
"""

import argparse

import os
from collections import defaultdict

from SemanticRoleChain_to_Echarts import generate_graph_option

from Type_Statistic_Bar_Generation import generate_type_sta_bar_option

from Gene_Related_Concept import generate_table_html

def read_html(html_frag_file: str):
    html_doc = open(html_frag_file).read()
    return html_doc

def main(chain_file: str, save_file: str, selected_gene_file: str,
         html_head: str='../data/html_fragments/html_head.txt'):

    """
    1. read html head file.
    2. get echarts option
    3. integrate html file.
    """

    graph_option = generate_graph_option(chain_file, only_option=True)

    type_sta_bar_option = generate_type_sta_bar_option(chain_file)

    table_html = generate_table_html(selected_gene_file)

    html_doc = read_html(html_head)

    html_doc = html_doc.replace('table_html_replace_here', table_html)

    html_doc += '<script type="text/javascript">\n'

    html_doc += type_sta_bar_option

    html_doc += '\n'

    html_doc += graph_option

    html_doc += """
    </script>
    </body>
    </html>
    """

    with open(save_file, 'w') as wf:
        wf.write(f'{html_doc}')

    print(f'{save_file} save done.')

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Web Integration.')
    parser.add_argument('-cf', dest='chain_file', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)

    parser.add_argument('-sg', dest='selected_gene_file', required=True)
    args = parser.parse_args()

    main(args.chain_file, args.save_file, args.selected_gene_file)
