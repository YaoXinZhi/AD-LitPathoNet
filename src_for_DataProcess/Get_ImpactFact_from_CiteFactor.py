# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 18/07/2022 10:20
@Author: XINZHI YAO
"""
import re

import requests

from bs4 import BeautifulSoup

def url_parser(url: str):
    html = requests.get(url)
    html = html.content.decode('utf-8')
    soup = BeautifulSoup(html, 'html.parser')
    return soup

def sub_parser(url: str):
    soup = url_parser(url)

    # impact_factor = 1
    circle_ele = soup.find('div', attrs={'class': 'circle-number'})
    if circle_ele:
        impact_factor = circle_ele.string
    else:
        print('Empty impact factor')
        print(f'{url}')
        impact_factor = 'None'
        # input()
    return impact_factor


def main_parser(url: str, save_file: str):
    soup = url_parser(url)

    body = soup.body

    ul = body.find('div', attrs={'class':'col-sm-12'}).find('ul')

    journal = ''
    sub_url = ''
    save_count = 0
    with open(save_file, 'w') as wf:
        wf.write(f'Journal\tISSN\tEISSN\tImpactFactor\n')
        for ele in ul.children:
            if ele.name == 'li':
                sub_url = ele.find('a').attrs['href']
                journal = ele.string
            if ele.name == 'span':
                if not journal or not sub_url:
                    print('empty journal.')
                    input()

                id_str = ele.string.split('/')
                if len(id_str) == 2:
                    issn, eissn = id_str
                else:
                    print('Wrong ID string.')
                    print(ele.string)
                    issn = id_str[0]
                    eissn = id_str[1]

                issn = re.findall(r'\d+-\d+', issn)
                if issn:
                    issn = issn[0]
                else:
                    issn = 'None'

                eissn = re.findall(r'\d+-\d+', eissn)
                if eissn:
                    eissn = eissn[0]
                else:
                    eissn = 'None'

                impact_factor = sub_parser(sub_url)

                wf.write(f'{journal}\t'
                         f'{issn}\t{eissn}\t'
                         f'{impact_factor}\n')

                journal = ''
                sub_url = ''

                save_count += 1
                if save_count % 50 ==0:
                    print(f'{save_count:,} Journals saved.')

    print(f'{save_file} save done.')

if __name__ == '__main__':

    # 20220718 CiteFactor
    CiteFactor_url = 'https://www.citefactor.org/journal-impact-factor-report-2021.html'

    IF_save_file = '../data/ImpactFactor.2022.tsv'

    main_parser(CiteFactor_url, IF_save_file)


