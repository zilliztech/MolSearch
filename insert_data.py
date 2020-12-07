import time
from milvus import *
import getopt
import os
import sys
import psycopg2
import time
from functools import reduce
import numpy as np
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import AllChem
from milvus import *

SERVER_ADDR = "192.168.1.85"
SERVER_PORT = 19530
VECTOR_DIMENSION = 2048


PG_HOST = "127.0.0.1"
PG_PORT = 5432
PG_USER = "postgres"
PG_PASSWORD = "postgres"
PG_DATABASE = "postgres"
PG_TABLE_NAME = "milvus"

MILVUS = Milvus(host=SERVER_ADDR, port=SERVER_PORT, timeout=100)


def connect_postgres_server():
    try:
        conn = psycopg2.connect(host=PG_HOST, port=PG_PORT, user=PG_USER, password=PG_PASSWORD, database=PG_DATABASE)
        print("connect the database!")
        return conn
    except:
        print("unable to connect to the database")


def create_pg_table(conn, cur):
    sql = "CREATE TABLE " + PG_TABLE_NAME + " (smiles text, ids bigint);"
    print(sql)
    try:
        cur.execute(sql)
        conn.commit()
        print("create postgres table!")
    except:
        print("can't create postgres table")
        # sys.exit()


def copy_data_to_pg(conn, cur, file_path):
    sql1 = "copy " + PG_TABLE_NAME + " from '" + file_path + "';"
    sql2 = "create index ids on " + PG_TABLE_NAME + "(ids);"
    print(sql1,sql2)
    try:
        cur.execute(sql1)
        conn.commit()
        cur.execute(sql2)
        conn.commit()
        print("copy data to pg sucessful!")
    except:
        print("faild  copy!")


def smiles_to_vec(smiles):
    mols = Chem.MolFromSmiles(smiles)
    # fp = AllChem.GetMorganFingerprintAsBitVect(mols, 2, VECTOR_DIMENSION)
    fp = Chem.RDKFingerprint(mols, fpSize=VECTOR_DIMENSION)
    hex_fp = DataStructs.BitVectToFPSText(fp)
    # print(hex_fp)
    vec = bytes.fromhex(hex_fp)
    return vec


def feature_extract(filepath):
    feats = []
    smiles = []
    ids = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.split()
            # print(parts)
            try:
                ids.append(int(parts[1]))
                smiles.append(parts[0])
                vec = smiles_to_vec(parts[0])
                feats.append(vec)
            except :
                print(line)
    return feats, smiles, ids



def do_load(file_path):
    vectors, names, ids = feature_extract(file_path)
    print("-----len of vectors:",len(vectors))
    connect_milvus_server()


    if table_name not in MILVUS.list_collections():
        collection_param = {
        "fields": [
            {"name": "embedding", "type": DataType.BINARY_VECTOR, "params": {"dim": VECTOR_DIMENSION}},
        ],
        "segment_row_limit": 800000,
        "auto_id": True
        }
        status = MILVUS.create_collection(param)
        print(status)


    ids_lens = 0
    while ids_lens<len(vectors) :
        try:
            status, ids = MILVUS.insert(collection_name=table_name, records=vectors[ids_lens:ids_lens+200000], ids = ids[ids_lens:ids_lens+200000])
        except:
            status, ids = MILVUS.insert(collection_name=table_name, records=vectors[ids_lens:len(vectors)], ids = ids[ids_lens:len(vectors)])
        ids_lens += 200000
        print(status, "ids:", len(ids))

    conn = connect_postgres_server()
    cur = conn.cursor()
    create_pg_table(conn, cur)
    copy_data_to_pg(conn, cur, file_path)
    cur.close()
    conn.close()


def main(argv):
    try:
        opts, args = getopt.getopt(
            sys.argv[1:],
            "f:h",
            ["help","file="],
        )
        # print(opts)
    except getopt.GetoptError:
        print("Usage:python insert_data.py --file <file_path> ")
        sys.exit(2)

    for opt_name, opt_value in opts:
        if opt_name in ("-f", "--file"):
            file_path = opt_value
            # connect_milvus_server()
            do_load(file_path)
        else:
            print("wrong parameter")
            sys.exit(2)

if __name__ == "__main__":
    main(sys.argv[1:])
