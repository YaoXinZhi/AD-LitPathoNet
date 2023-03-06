# AD-LitPathoNet
**AD-LitPathoNet: A Resource of Pathology Network with Rich Literature Evidence for Alzheimer’s Disease**

Pathology Network with Rich Literature Evidence for Alzheimer's disease (AD-LitPathoNet) is a text-mining based database, which is aim to facilitate the systematic knowledge curation and pathology mechanism investigation for AD.

AD-LitPathoNet keeps quarterly updates, and the latest data will updated to the database, which you can download directly without running any code.

## System Requirements  
The project code is mainly written in python 3.7 and depends on several APIs such as E-Utils, PubTator API, OGER, Phenotype. The complete pipeline dependencies can be found in requirement.txt.

The code was tested on centos 7 and can also be run on other linux kernel systems or windows. Most of the pipeline code supports multi-threading, which means more cpu cores to help update data faster, and a smooth network is needed to access these APIs.

AD-LitPathoNet relies in part on several other projects for data annotation, including the AGAC Biomedical Entity NER (https://github.com/YaoXinZhi/BERT-CRF-for-BioNLP-OST2019-AGAC-Task1) and the AGAC Thesis Role Identification ( https://github.com/YaoXinZhi/BERT-for-BioNLP-OST2019-AGAC-Task2).

