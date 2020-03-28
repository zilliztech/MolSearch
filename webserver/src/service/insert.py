import logging
from common.const import default_cache_dir
from common.config import PG_TABLE
from indexer.index import milvus_client, create_table, insert_vectors, delete_table, search_vectors, create_index
from diskcache import Cache
from encoder.encode import smiles_to_vec
import psycopg2


PG_HOST = "localhost"
PG_PORT = 5432
PG_USER = "zilliz"
PG_PASSWORD = "zilliz"
PG_DATABASE = "milvus"


def connect_postgres_server():
    try:
        conn = psycopg2.connect(host=PG_HOST, port=PG_PORT, user=PG_USER, password=PG_PASSWORD, database=PG_DATABASE)
        print("connect the database!")
        return conn
    except:
        print("unable to connect to the database")


def search_loc_in_pg(cur, table_name=PG_TABLE):
    try:
        sql = "select smiles from " + table_name+ " where ids = '" + str(ids) + "';"
        cur.execute(sql)
        rows = cur.fetchall()
        return str(rows[0][0])
    except:
        print("search faild!")


def main():
    conn = connect_postgres_server()
    cur = conn.cursor()
    insert_data_to_pg(cur)


if __name__ == '__main__':
    main()