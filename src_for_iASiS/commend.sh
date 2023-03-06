# 统计链条文件中概念 用于iASiS数据检索
Concept_statistic.py
Concept_to_umls.py
GO_HPO_MESH_UMLS_mapping.py
# 开启Neo4j 下载iASiS数据
iasis_search.py
Bulk_iASiS_query.py
# iASiS数据转换为PNRLE标准化体系
iASiS_to_PNRLE.py
# iASiS数据根据55中SemanticRelation进行筛选
iASiS_filter.py
# 利用networkx包结合两套数据进行通路搜索
Pathway_seeker.py
# pathway中每个实体数量统计及存在实体名字相同 id不同实体的统计
pathway_info_statistic.py
# 基于规则和先验知识对通路进行质控，打分
pathway_filter.py
# 整理iASiS的边和证据信息为Django数据库格式 用于后端搜索
iASiS_and_pathway_to_DB.py
