import getopt
import os
import sys
import time
from functools import reduce

import numpy as np
from milvus import *

SERVER_ADDR = "192.168.1.85"
SERVER_PORT = 19530

FILE_NPY_PATH = '/data/pubchem/out_data/out_npy'
FILE_IDS = '/data/pubchem/out_data/out_ids'


milvus = Milvus()

# is_uint8 = False
# if_normaliz = False


def normaliz_data(vec_list):
    for i in range(len(vec_list)):
        vec = vec_list[i]
        square_sum = reduce(lambda x, y: x + y, map(lambda x: x * x, vec))
        sqrt_square_sum = np.sqrt(square_sum)
        coef = 1 / sqrt_square_sum
        vec = list(map(lambda x: x * coef, vec))
        vec_list[i] = vec
    return vec_list


def load_npy_data(filename):
    filename = FILE_NPY_PATH + "/" + filename
    print("npy_file:",filename)
    data = np.load(filename) 
    if is_uint8:
        data = (data+0.5)/255
    if if_normaliz:
        data = normaliz_data(data)
    data = data.tolist()
    return data


def load_uint8_data(filename):
    filename = FILE_UINT8_PATH + "/" + filename
    print("uint8_file:",filename)
    data = []
    for line in open(filename, 'r'):
        line = line.strip('\n')
        data_uint8 = line.split(' ')
        data_uint8 = list(map(int, data_uint8))
        data_bytes = bytes(data_uint8)
        # print(data_uint8,data_bytes,'\n')
        data.append(data_bytes)
    return data


def load_csv_data(filename):
    import pandas as pd
    filename = FILE_CSV_PATH + "/" + filename
    # print(filename)
    data = pd.read_csv(filename, header=None)
    data = np.array(data)
    if is_uint8:
        data = (data+0.5)/255
    if if_normaliz:
        data = normaliz_data(data)
    data = data.tolist()
    # data = normaliz_data(data)
    return data


def load_fvecs_data(fname, base_len, idx):
    begin_num = base_len * idx
    # print(fname, ": ", begin_num )
    x = np.memmap(fname, dtype='uint8', mode='r')
    d = x[:4].view('int32')[0]
    data = x.view('float32').reshape(-1, d + 1)[begin_num:(begin_num + base_len), 1:]
    if is_uint8:
        data = (data+0.5)/255
    if if_normaliz:
        data = normaliz_data(data)
    data = data.tolist()
    return data


def load_bvecs_data(fname, base_len, idx):
    begin_num = base_len * idx
    # print(fname, ": ", begin_num)
    x = np.memmap(fname, dtype='uint8', mode='r')
    d = x[:4].view('int32')[0]
    data = x.reshape(-1, d + 4)[begin_num:(begin_num + base_len), 4:]
    if is_uint8:
        data = (data+0.5)/255
    if if_normaliz:
        data = normaliz_data(data)
    data = data.tolist()
    return data


def handle_status(status):
    if status.code != Status.SUCCESS:
        print(status)
        sys.exit(2)


def connect_milvus_server():
    print("connect to milvus")
    status = milvus.connect(host=SERVER_ADDR, port=SERVER_PORT, timeout=1000 * 1000 * 20)
    handle_status(status=status)
    return status


"""
def npy_to_milvus(MILVUS_TABLE):
    filenames = os.listdir(FILE_NPY_PATH)
    filenames.sort()

    filenames_ids = os.listdir(FILE_IDS)
    filenames_ids.sort()
    count = 0
    for filename in filenames:
        ids_vec = load_ids(filenames_ids[count])
        vectors = load_hex(filename)
        # vectors = load_npy_data(filename)
        ids_lens = 0
        print(len(ids_vec),len(vectors),'\n')
        time_add_start = time.time()
        # status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors, ids=ids_vec)
        status, ids = milvus.insert(collection_name=MILVUS_TABLE, records=vectors, ids=ids_vec)
        print(status)
        time_add_end = time.time()
        print(filename, " insert milvus time-------: ", time_add_end - time_add_start)
        count += 1
"""

def npy_to_milvus(MILVUS_TABLE):
    filenames = os.listdir(FILE_NPY_PATH)
    filenames.sort()

    # filenames_ids = os.listdir(FILE_IDS)
    # filenames_ids.sort()
    count = 0
    for filename in filenames:
        if count>10:
            break
        # ids_vec = load_ids(filenames_ids[count])
        ids_vec = []
        vectors = load_hex(filename)
        for i in range(len(vectors)):
            location = '8' + '%05d'%count  + '%07d'%i
            ids_vec.append(int(location))
        # vectors = load_npy_data(filename)
        ids_lens = 0
        print(len(ids_vec),len(vectors),'\n')
        while ids_lens<len(vectors) :
            time_add_start = time.time()
            try:
                # status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors[ids_lens:ids_lens+100000], ids=ids_vec[ids_lens:ids_lens+100000])
                status, ids = milvus.insert(collection_name=MILVUS_TABLE, records=vectors[ids_lens:ids_lens+200000], ids=ids_vec[ids_lens:ids_lens+200000])
                print(status)
            except:
                # status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors[ids_lens:len(vectors)], ids=ids_vec[ids_lens:len(vectors)])
                status, ids = milvus.insert(collection_name=MILVUS_TABLE, records=vectors[ids_lens:len(vectors)], ids=ids_vec[ids_lens:len(vectors)])
                print(status)
            time_add_end = time.time()
            print("ids:",len(ids),ids_vec[0])
            print(filename, " insert milvus time: ", time_add_end - time_add_start)
            ids_lens += 200000
        count += 1



def uint8_to_milvus(MILVUS_TABLE):
    filenames = os.listdir(FILE_UINT8_PATH)
    filenames.sort()

    filenames_ids = os.listdir(FILE_IDS)
    filenames_ids.sort()
    # file_index = 0
    count = 0
    for filename in filenames:
        ids_vec = load_ids(filenames_ids[count])
        vectors = load_uint8_data(filename)
        ids_lens = 0
        print("ids:", len(ids_vec), "vectors:", len(vectors))
        while ids_lens<len(vectors) :
            time_add_start = time.time()
            # print(vectors_ids)
            try:
                status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors[ids_lens:ids_lens+100000], ids=ids_vec[ids_lens:ids_lens+100000])
            except:
                status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors[ids_lens:len(vectors)], ids=ids_vec[ids_lens:len(vectors)])
            print(status)
            time_add_end = time.time()
            # print("ids:",len(ids),ids_vec[0])
            print("insert milvus time: ", time_add_end - time_add_start, '\n')
            ids_lens += 100000
            # file_index = file_index + 1
        count += 1


def csv_to_milvus(MILVUS_TABLE):
    filenames = os.listdir(FILE_CSV_PATH)
    filenames.sort()
    file_index = 0
    for filename in filenames:
        vectors = load_csv_data(filename)
        vectors_ids = []
        for i in range(len(vectors)):
            location = '8' + '%04d'%file_index  + '%06d'%i
            vectors_ids.append(int(location))
        time_add_start = time.time()
        status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors, ids=vectors_ids)
        time_add_end = time.time()
        print(filename, " insert time: ", time_add_end - time_add_start)
        file_index = file_index + 1


def load_ids(file):
    file = FILE_IDS + '/' + file
    print("ids_file:",file)
    ids = []
    for line in open(file, 'r'):
        data = line.strip('\n')
        #ids.append(int(data[4:]))
        ids.append(int(data))
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


def main(argv):
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "ncfbt:u",
            ["npy", "uint8", "csv", "fvecs", "bvecs","table="],
        )
        # print(opts)
    except getopt.GetoptError:
        print("Usage: python milvus_load.py -t <table> -n")
        sys.exit(2)

    for opt_name, opt_value in opts:
        if opt_name in ("-t", "--table"):
            MILVUS_TABLE = opt_value
            PG_TABLE_NAME = opt_value
        elif opt_name in ("-n", "--npy"):
            connect_milvus_server()
            npy_to_milvus(MILVUS_TABLE)
        elif opt_name in ("-u", "--u8"):
            connect_milvus_server()
            uint8_to_milvus(MILVUS_TABLE)
        elif opt_name in ("-c", "--csv"):
            connect_milvus_server()
            csv_to_milvus(MILVUS_TABLE)
        elif opt_name in ("-f", "--fvecs"):
            connect_milvus_server()
            count = 0
            while count < (FVECS_VEC_NUM // FVECS_BASE_LEN):
                vectors = load_fvecs_data(FILE_FVECS_PATH, FVECS_BASE_LEN, count)
                print(count*FVECS_BASE_LEN, " ", (count+1)*FVECS_BASE_LEN)
                vectors_ids = [id for id in range(count*FVECS_BASE_LEN,(count+1)*FVECS_BASE_LEN)]                
                time_add_start = time.time()
                status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors, ids=vectors_ids)
                time_add_end = time.time()
                print(count, " insert to milvus time: ", time_add_end - time_add_start)
                count = count + 1
                
        elif opt_name in ("-b", "--bvecs"):
            connect_milvus_server()
            count = 0
            while count < (FVECS_VEC_NUM // FVECS_BASE_LEN):
                vectors = load_bvecs_data(FILE_BVECS_PATH, FVECS_BASE_LEN, count)
                print(count*FVECS_BASE_LEN, " ", (count+1)*FVECS_BASE_LEN)
                vectors_ids = [id for id in range(count*FVECS_BASE_LEN,(count+1)*FVECS_BASE_LEN)]                
                time_add_start = time.time()
                status, ids = milvus.add_vectors(table_name=MILVUS_TABLE, records=vectors, ids=vectors_ids)
                time_add_end = time.time()
                print(count, " insert to milvus time: ", time_add_end - time_add_start)
                count = count + 1
        else:
            print("wrong parameter")
            sys.exit(2)

if __name__ == "__main__":
    main(sys.argv[1:])


