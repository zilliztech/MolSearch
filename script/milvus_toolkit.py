import time
import os
import getopt
import sys
import datetime
import numpy as np
from milvus import *


MILVUS = Milvus()
SERVER_ADDR = "192.168.1.85"
SERVER_PORT = 19530

index_file_size = 2048
# metric_type = MetricType.SUBSTRUCTURE
metric_type = MetricType.SUPERSTRUCTURE
# metric_type = MetricType.JACCARD
nlist=2048

NL_FOLDER_NAME = '/data/source_data'
NQ_FOLDER_NAME = 'nq_npy_2048'
PE_FOLDER_NAME = 'performance'

nq_scope = [1]
topk_scope = []
for i in range(50):
    topk_scope.append(100)


# the status of milvus
def handle_status(status):
    if status.code != Status.SUCCESS:
        print(status)
        sys.exit(2)


# connect to the milvus server
def connect_server():
    print("connect to milvus")
    status = MILVUS.connect(host=SERVER_ADDR, port=SERVER_PORT, timeout=100)
    handle_status(status=status)
    return status


# -c/create the table with milvus
def create_table(table_name, dim):
    param = {'collection_name': table_name, 'dimension': dim, 'index_file_size':index_file_size, 'metric_type':metric_type}
    print("create table: ", table_name, " dimension:", dim, 'metric_type:', metric_type)
    return MILVUS.create_collection(param)


def table_show():
    print(MILVUS.show_collections())


def describe_table(table_name):
    print(MILVUS.describe_collection(table_name)[1])


def delete_table(table_name):
    print("delete table:", table_name)
    MILVUS.drop_collection(table_name)


def build_table(table_name,index_type):
    # index_param = {'index_type': index_type, 'nlist': nlist}
    # status = MILVUS.create_index(table_name, index_param)
    index_param = {'nlist': nlist}
    status = MILVUS.create_index(table_name, index_type, index_param)
    print(status)


def show_server_version():
    print(MILVUS.server_version()[1])


def show_client_version():
    print(MILVUS.client_version())


def table_rows(table_name):
    print(table_name, 'has', MILVUS.count_collection(table_name)[1], 'rows')


def load_nq_vec(nq):
    vectors = []
    length = 0
    filenames = os.listdir(NQ_FOLDER_NAME)
    filenames.sort()
    for filename in filenames:
        vec_list = load_vec_list(NQ_FOLDER_NAME + '/' + filename)
        length += len(vec_list)
        if length > nq:
            num = nq % len(vec_list)
            vec_list = vec_list[0:num]
        vectors += vec_list
        if len(vectors) == nq:
            return vectors


def load_mols_vec(nq):
    vectors = []
    length = 0
    filenames = os.listdir(NQ_FOLDER_NAME)
    filenames.sort()
    for filename in filenames:
        # vec_list = load_uint8_vec(NQ_FOLDER_NAME + '/' + filename)
        vec_list = load_hex_vec(NQ_FOLDER_NAME + '/' + filename)
        length += len(vec_list)
        if length > nq:
            num = nq % len(vec_list)
            vec_list = vec_list[0:num]
        vectors += vec_list
        if len(vectors) == nq:
            return vectors


def load_uint8_vec(filename):
    data = []
    for line in open(filename, 'r'):
        line = line.strip('\n')
        data_uint8 = line.split(' ')
        data_uint8 = list(map(int, data_uint8))
        data_bytes = bytes(data_uint8)
        # print(data_uint8,data_bytes,'\n')
        data.append(data_bytes)
    return data


def load_hex_vec(filename):
    data = np.load(filename)
    data = data.tolist()
    vectors = []
    for d in data:
        vectors.append(bytes.fromhex(d))
    return vectors


def load_vec_list(file_name, num=0):
    if IS_CSV:
        import pandas as pd
        data = pd.read_csv(file_name, header=None)
        data = np.array(data)
    else:
        data = np.load(file_name)
    if IS_UINT8:
        data = (data + 0.5) / 255
    vec_list = data.tolist()
    return vec_list


def search_mols_list(table_name,np):
    # random1 = datetime.datetime.now().strftime("%m%d%H%M")
    if not os.path.exists(PE_FOLDER_NAME):
        os.mkdir(PE_FOLDER_NAME)
    filename = PE_FOLDER_NAME + '/' + table_name + '_' + str(np) +".txt" 
    file = open(filename, "w+")
    file.write('nq,topk,total_time,avg_time' + '\n')
    for nq in nq_scope:
        time_start = time.time()
        query_list = load_mols_vec(nq)
        time_end = time.time()
        print("load query:", len(query_list), "time_load = ", time_end - time_start)
        for k in topk_scope:
            time_start = time.time()
            # MILVUS.search_vectors(table_name=table_name, query_records=query_list, top_k=k, nprobe=np)
            MILVUS.search(collection_name=table_name, query_records=query_list, top_k=k, params={"nprobe": np})
            time_end = time.time()
            time_cost = time_end - time_start
            line = str(nq) + ',' + str(k) + ',' + str(round(time_cost, 4)) + ',' + str(round(time_cost / nq, 4)) + '\n'
            file.write(line)
            print(nq, k, time_cost)
        file.write('\n')
    file.close()
    print("search_vec_list done !")


def search_vec_list(table_name,np):
    # random1 = datetime.datetime.now().strftime("%m%d%H%M")
    if not os.path.exists(PE_FOLDER_NAME):
        os.mkdir(PE_FOLDER_NAME)
    filename = PE_FOLDER_NAME + '/' + table_name + '_' + str(np) + PE_FILE_NAME
    file = open(filename, "w+")
    file.write('nq,topk,total_time,avg_time' + '\n')
    for nq in nq_scope:
        time_start = time.time()
        query_list = load_nq_vec(nq)
        time_end = time.time()
        print("load query:", len(query_list), "time_load = ", time_end - time_start)
        for k in topk_scope:
            time_start = time.time()
            MILVUS.search_vectors(table_name=table_name, query_records=query_list, top_k=k, nprobe=np)
            time_end = time.time()
            time_cost = time_end - time_start
            line = str(nq) + ',' + str(k) + ',' + str(round(time_cost, 4)) + ',' + str(round(time_cost / nq, 4)) + '\n'
            file.write(line)
            print(nq, k, time_cost)
        file.write('\n')
    file.close()
    print("search_vec_list done !")


def is_normalized():
    filenames = os.listdir(NL_FOLDER_NAME)
    filenames.sort()
    vetors = load_vec_list(NL_FOLDER_NAME+'/'+filenames[0])
    for i in range(10):
        sqrt_sum = np.sum(np.power(vetors[i], 2))
        print(sqrt_sum)

def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "chdsn",
            ["help", "table=", "dim=", "index=", "nq=", "show", "describe", "delete", "build", "drop_index", "server_version",
             "client_version", "rows", "normal", "nprobe=", "has", "desc_index"]
        )
    except getopt.GetoptError:
        print("Usage: python milvus_toolkit.py -t <table> -d <dimension> -c")
        sys.exit(2)

    for opt_name, opt_value in opts:
        if opt_name in ("--help"):
            print("python milvus_toolkit.py -t <table> --build")
            sys.exit()
        elif opt_name == "--table":
            table_name = opt_value
        elif opt_name == "--dim":
            dim = int(opt_value)
        elif opt_name == "--index":
            index_type = opt_value
            if index_type == 'flat':
                it = IndexType.FLAT
            elif index_type == 'ivf':
                it = IndexType.IVFLAT
            elif index_type == 'sq8':
                it = IndexType.IVF_SQ8
        elif opt_name == "-c":
            connect_server()
            status = create_table(table_name, dim)
            print(status)
        elif opt_name == "--show":
            connect_server()
            table_show()
        elif opt_name in ("-n", "--normal"):
            is_normalized()
        elif opt_name == "--describe":
            connect_server()
            describe_table(table_name)
        elif opt_name == "--has":
            connect_server()
            print(MILVUS.has_table(table_name=table_name))
        elif opt_name in ("-d", "--delete"):
            connect_server()
            delete_table(table_name)
        elif opt_name == "--build":
            connect_server()
            print(table_name," ",it)
            time1 = time.time()
            build_table(table_name,it)
            time2 = time.time()
            print("total cost time: ", time2-time1)
        elif opt_name == "--drop_index":
            connect_server()
            print(MILVUS.drop_index(table_name))
        elif opt_name == "--desc_index":
            connect_server()
            print(MILVUS.describe_index(table_name))
        elif opt_name == "--server_version":
            connect_server()
            show_server_version()
        elif opt_name == "--client_version":
            connect_server()
            show_client_version()
        elif opt_name == "--rows":  # test.py --table <tablename> --rows
            connect_server()
            table_rows(table_name)
        elif opt_name == "--nprobe":
            np = int(opt_value)
        elif opt_name == "-s":  # test.py --table <tablename> --nprobe <np> -s
            connect_server()
            # search_vec_list(table_name,np)
            search_mols_list(table_name,np)
            sys.exit()


if __name__ == '__main__':
    main()
