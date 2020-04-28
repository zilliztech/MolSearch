# README

## insert_data.py说明

### 参数说明
| 参数             | 描述                        | 默认设置     |
| ---------------- | --------------------------- | ------------ |
| SERVER_ADDR      | milvus server 的 IP 地址    | 192.168.1.85 |
| SERVER_PORT      | milvus server 端口号        | 19530        |
| VECTOR_DIMENSION | 向量维度                    | 2048         |
| table_names      | 在milvus中的表              | [...]        |
| metric_types     | 对应milvus表中不同的metrics | [...]        |
| PG_HOST          | postgres 的 IP 地址         | 192.168.1.85 |
| PG_PORT          | postgres 的端口号           | 5430         |
| PG_USER          | postgres 的用户名           | postgres     |
| PG_PASSWORD      | postgres 的密码             | postgres     |
| PG_DATABASE      | postgres 的数据库名         | postgres     |
| PG_TABLE_NAME    | postgres 的表名             | milvus       |

### 使用说明

```bash
$ python insert_data -f '<file_path>'
# 执行此命令讲文件中的smiles转为fingerprint并导入milvus,并在pg中建立对应关系表
```





## thread_gen_smiles_fp.py说明

### 参数说明

| 参数        | 描述                                       | 默认设置                 |
| ----------- | ------------------------------------------ | ------------------------ |
| file_length | 每个文件保存 fingerprint 的数量            | 200000                   |
| vec_dim     | 生成 fingerprint 的维度/位数               | 2048                     |
| FILE_PATH   | 保存.smi文件的目录                         | /data/workspace/test     |
| OUT         | 生成结果的目录                             | /data/workspace/out_test |
| OUT_SMILES  | 转 fingerprint 成功的 smiles 文件路径      | out_smiles               |
| OUT_IDS     | 与 smiles 对应的 ID 文件路径               | out_ids                  |
| OUT_NPY     | 与 smiles 对应 fingerprint 的 npy 文件路径 | out_npy                  |

### 使用说明

```bash
$ python thread_gen_smiles_fp.py
# 执行此命令使用多线程将 FILE_PATH 下的 smiles 转化为 RDKit fingerprint，程序将建立 OUT 等目录存储结果。
```



## milvus_toolkit.py说明

### 参数说明：

| 参数            | 描述                            | 默认设置           |
| --------------- | ------------------------------- | ------------------ |
| SERVER_ADDR     | milvus server 的 IP 地址        | 192.168.1.85       |
| SERVER_PORT     | milvus server 端口号            | 19530              |
| index_file_size | milvus 建立索引时的文件大小阈值 | 512                |
| metric_type     | milvus search的类型             | MetricType.JACCARD |
| nlist           | milvus search时所分的桶数量     | 2048               |

### 使用说明：

```bash
$python milvus_toolkit.py --table <table_name> --dim <dim_num> -c
# 执行-c，在milvus中建表
# -t或者--table表示表名
# -d或者--dim 表示维度/位数

$ python3 milvus_toolkit.py --show
# 执行--show，显示milvus中所有表的表名

$ python milvus_toolkit.py --table <table_name> --index ivf --build
# 执行--build，给表建立IVFLAT索引
# -t或者--table表示表名
# --index 表示索引类型
```



## milvus_load.py说明

### 参数说明

| 参数          | 描述                                | 默认设置                         |
| ------------- | ----------------------------------- | -------------------------------- |
| SERVER_ADDR   | milvus server 链接地址              | 192.168.1.85                     |
| SERVER_PORT   | milvus server端口号                 | 19530                            |
| FILE_NPY_PATH | 导入数据时的 npy 向量所在文件夹路径 | /data/workspace/out_test/out_npy |
| FILE_IDS      | 导入数据时的 ids 所在文件夹路径     | /data/workspace/out_test/out_ids |

### 使用说明

```bash
$ python milvus_load.py --table <table_name> -n
# 执行-n，将存储格式为npy的向量导入milvvus
# -t或者--table表示表名
```



## milvus_search.py说明

### 参数说明

| 参数                      | 描述                       | 默认设置       |
| ------------------------- | -------------------------- | -------------- |
| SERVER_ADDR               | Milvus 的IP设置            | 192.168.1.85 |
| SERVER_PORT               | Milvus 的端口设置          | 19530          |
| NQ_FOLDER_NAME            | 查询向量集的路径           | nq_npy |
| SE_FOLDER_NAME            | 查询结果保存的路径         | search       |
| SE_FILE_NAME              | 查询结果保存的文件名       | _output.txt  |
| BASE_FOLDER_NAME (ignore) | 源向量数据集的路径         | /data/milvus |
| TOFILE(ignore)            | 是否存储查询后的文件信息   | True           |
| **GT_NQ**                 | **ground truth中的nq数值** | **100** |
| NPROBE                    | Milvus参数nprobe           | 64          |

### 使用说明

```bash
$ python milvus_search.py -table <tablename> -q <nq> -k <topk> -n <nprobe> -s

# 执行-s实现Milvus的向量查询，并将结果写入SEARCH_FOLDER_NAME目录下的table_name_output.txt中，该文件有随机数，查询结果ids和查询结果distance三列
# -t或者--table表示需要查询的表名
# -q或者--nq表示在查询集中随机选取的查询向量个数，该参数可选，若没有-q表示查询向量为查询集中的全部数据
# -k或者--topk表示查询每个向量的前k个相似的向量
# -n或者--nprobe表示milvus参数NPROBE
```



## get_results_smiles.py说明

### 参数说明

| 参数            | 描述                          | 默认设置                            |
| --------------- | ----------------------------- | ----------------------------------- |
| SE_FOLDER_NAME  | 查询结果保存的路径            | search                              |
| SE_FILE_NAME    | 查询结果保存的文件名          | _output.txt                         |
| CM_FOLDER_NAME  | 生成对应smiles的路径          | compare                             |
| CM_GET_LOC_NAME | 生成对应smiles的文件名        | _compare.txt                        |
| FILE_SMILES     | 与 ids 对应的 smiles 文件路径 | /data/workspace/out_test/out_smiles |

### 使用说明

```bash
$ python get_results_smiles.py --table=<table_name> -n <nprobe> -g
# 执行-g生成与milvus结果对应的smiles
# -t或者--table表示表名
# -n或者--nprobe表示milvus参数NPROBE
```
