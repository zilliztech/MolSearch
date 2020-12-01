import logging
import time
from common.config import DEFULT_TABLE
from diskcache import Cache
from indexer.index import milvus_client, count_table


def do_count(table_name):
    if not table_name:
        table_name = DEFULT_TABLE
    try:
        index_client = milvus_client()
        print("get table rows:",table_name)
        num = count_collection(index_client, table_name)
        return num
    except Exception as e:
        logging.error(e)
        return "Error with {}".format(e)