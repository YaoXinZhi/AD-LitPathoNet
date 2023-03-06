# -*- coding:utf-8 -*-
# ! usr/bin/env python3
"""
Created on 14/06/2022 9:33
@Author: XINZHI YAO
"""

import argparse

from py2neo import Graph, NodeMatcher, RelationshipMatcher


def search(keyword: str, save_file: str):
    print(f'Query keyword: {keyword}.')

    graph = Graph("http://127.0.0.1:7474", auth=("neo4j", "happy19970308"))
    # node_label_set = graph.schema.node_labels
    # relation_set = {relation for relation in graph.schema.relationship_types if not relation.endswith('__SPEC__')}
    node_matcher = NodeMatcher(graph)
    relation_matcher = RelationshipMatcher(graph)

    query_node = node_matcher.match('Entity').where(label=keyword).first()

    relationships = list(relation_matcher.match([ query_node ], r_type=None))

    print(f'Query count: relation-{len(relationships)}')

    if len(relationships) != 0:
        with open(save_file, 'w') as wf:
            wf.write('Relation\tSource Name\tSource\tTarget Name\tTarget\n')
            for relation in relationships:
                start = relation.start_node
                start_label = start[ 'label' ]
                end = relation.end_node
                end_label = end[ 'label' ]

                rel = type(relation).__name__

                wf.write(f'{rel}\t'
                         f'{start_label}\t{str(start)}\t'
                         f'{end_label}\t{str(end)}\n')
        print(f'{save_file} save done.')
    else:
        print('The query is empty.')


def graph_search(keyword: str, save_file: str, save_mention_in: bool):
    print(f'Query keyword: {keyword}.')

    graph = Graph("http://127.0.0.1:7474", auth=("neo4j", "happy19970308"))

    # query_data = graph.run(f"MATCH (a)-[r]-(b {label:'{keyword}'}) return *")
    # query_data = graph.run(f"MATCH (a)-[r]-(b) return *")
    # query_data = graph.run("MATCH (a)-[r]-(b {label:'Apoptotic'}) return *")
    commend = 'MATCH (a)-[r]-(b {id:"' + keyword + '"}) return *'
    #print(commend)
    query_data = graph.run(commend)

    with open(save_file, 'w') as wf:
        wf.write('Relation Type\tRelation\t'
                 'Source label\tSource type\tSource SemanticType\tSource\t'
                 'Target label\tTarget type\tTarget SemanticType\tTarget\n')
        for node in query_data.data():

            source = node[ 'a' ]
            target = node[ 'b' ]
            relation = node[ 'r' ]

            source_type = str(source.labels)[ 1: ]
            source_label = source[ 'label' ]
            source_sem_types = source[ 'sem_types' ]

            target_type = str(target.labels)[ 1: ]
            target_label = target[ 'label' ]
            target_sem_types = target[ 'sem_types' ]


            relation_type = list(relation.types())[0]
            # fixme: save complete relationship

            if not save_mention_in and relation_type == 'MENTIONED_IN':
                continue

            wf.write(f'{relation_type}\t{relation}\t'
                     f'{source_label}\t{source_type}\t{source_sem_types}\t{str(source)}\t'
                     f'{target_label}\t{target_type}\t{target_sem_types}\t{str(target)}\n')
    print(f'{save_file} save done.')


def main():
    parser = argparse.ArgumentParser('iasis search.')
    parser.add_argument('-kw', dest='keyword', required=True)
    parser.add_argument('-sf', dest='save_file', required=True)

    parser.add_argument('-sm', dest='save_mention_in', action='store_false',
                        default=True)
    args = parser.parse_args()

    # search(args.keyword, args.save_file)
    graph_search(args.keyword, args.save_file, args.save_mention_in)


if __name__ == '__main__':
    main()
