import os

MILVUS_HOST = os.getenv("MILVUS_HOST", "192.168.1.58")
MILVUS_PORT = os.getenv("MILVUS_PORT", 19520)
VECTOR_DIMENSION = os.getenv("VECTOR_DIMENSION", 512)
NUM = os.getenv("TOPK_NUM", 100)
SIM_TABLE = os.getenv("SIM_TABLE", "similarity_table")
SUB_TABLE = os.getenv("SUB_TABLE", "substructure_table")
SUPER_TABLE = os.getenv("SUPER_TABLE", "superstructure_table")
PG_TABLE = os.getenv("PG_TABLE", "milvus_table")
