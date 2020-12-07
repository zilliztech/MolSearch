import getopt
import os
import sys
import time
from functools import reduce

import numpy as np
from milvus import *

SERVER_ADDR = "127.0.0.1"
SERVER_PORT = 19530
MILVUS_TABLE = "milvus_mols" 

# FILE_NPY_PATH = '/mnt/test/apptec_milvus/molview_demo/out_data_2048/out_npy'
# FILE_IDS = '/mnt/test/apptec_milvus/molview_demo/out_data_2048/out_ids'
FILE_NPY_PATH = '/workspace/test_data/out_npy'
FILE_IDS = '/workspace/test_data/out_ids'

milvus = Milvus(host=SERVER_ADDR, port=SERVER_PORT)


def npy_to_milvus(table_name = MILVUS_TABLE):
    filenames = os.listdir(FILE_NPY_PATH)
    filenames.sort()

    filenames_ids = os.listdir(FILE_IDS)
    filenames_ids.sort()

    for filename in filenames:
        ids_vec = load_ids(filenames_ids[count])
        vectors = load_hex(filename)
        print(len(ids_vec),len(vectors),'\n')
        
        hybrid_entities = [{"name": "embedding", "values": vectors, "type": DataType.BINARY_VECTOR}]
        ids = milvus.insert(table_name, hybrid_entities, ids_vec)
            
        time_add_end = time.time()
        print("ids:",len(ids))
        print(filename, " insert milvus time: ", time_add_end - time_add_start)


def load_ids(file):
    file = FILE_IDS + '/' + file
    print("ids_file:",file)
    ids = []
    for line in open(file, 'r'):
        data = line.strip('\n')
        ids.append(int(data[4:]))
        # ids.append(int(data))
    return ids

def load_hex(file):
    file = FILE_NPY_PATH + '/' + file
    print("hex_file:",file)
    data = np.load(file)
    data = data.tolist()
    vectors = []
    for d in data:
        vectors.append(bytes.fromhex(d))
    return vectors


def main():
    npy_to_milvus()

if __name__ == "__main__":
    main()