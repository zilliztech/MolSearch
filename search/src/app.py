import os
import os.path as path
import logging
from common.config import SIM_TABLE,SUB_TABLE,SUPER_TABLE,NUM
from common.const import UPLOAD_PATH
from service.search import do_search
from service.count import do_count
from indexer.index import milvus_client, create_table, insert_vectors, delete_table, search_vectors, create_index
from flask_cors import CORS
from flask import Flask, request, send_file, jsonify
from flask_restful import reqparse
from werkzeug.utils import secure_filename
import numpy as np
from numpy import linalg as LA
import time


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_PATH
app.config['JSON_SORT_KEYS'] = False
CORS(app)

model = None

@app.route('/api/v1/count')
def do_count_api():
    args = reqparse.RequestParser(). \
        add_argument('Table', type=str). \
        parse_args()
    table_name = args['Table']
    rows = do_count(table_name)
    return "{}".format(rows)


@app.route('/api/v1/search')
def do_search_api():
    args = reqparse.RequestParser(). \
        add_argument("Molecular", type=str). \
        add_argument("Type", type=str, default='similarity'). \
        add_argument("Cid", type=int). \
        parse_args()

    molecular_name = args['Molecular']
    search_type = args['Type']
    cid_num = args['Cid']

    if search_type == 'similarity':
        table_name = SIM_TABLE
    elif search_type == 'substructure':
        table_name = SUB_TABLE
    elif search_type == 'superstructure':
        table_name = SUPER_TABLE
    top_k = NUM

    if cid_num:
        from pubchempy import get_compounds, Compound
        comp = Compound.from_cid(cid_num)
        molecular_name = comp.isomeric_smiles

    if not molecular_name:
        return "no molecular"
    if molecular_name:
        try:
            res_smi = do_search(table_name, molecular_name, top_k)
        except:
            return "There has no results, please input the correct molecular and ensure the table has data."
        re = {}
        re["Smiles"] = res_smi
        return jsonify(re), 200
    return "not found", 400


if __name__ == "__main__":
    app.run(host="0.0.0.0")
