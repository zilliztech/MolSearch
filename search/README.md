# Search molecular

## Api
### /api/v1/search 
#### methods
    POST
#### PARAM
||||
|-|-|-|
|Table|str|milvus table,defult milvus|
|Num|int|top k|
|Molecular|str|COc1ccc(cc1)SCCC(=O)NCCNS(=O)(=O)c1cccc(c1)Cl|

### /api/v1/count
#### methods
	POST
#### PARAM
||||
|-|-|-|
|Table|str|milvus table,defult milvus|



# Requirements

    milvus container



# Env

|||
|-|-|
|MILVUS_HOST |milvus container host|
|MILVUS_PORT |milvus container port|
|VECTOR_DIMENSION |default vector dimension number|
|DATA_PATH |mols data path|
|SIM_TABLE |similarity_table|
|SUB_TABLE |superstructure_table|
|SUPER_TABLE |substructure_table|



## how to use

    python3 src/app.py
