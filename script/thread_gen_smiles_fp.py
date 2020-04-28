import threading
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import ProcessPoolExecutor
import os
import time
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import Draw
import math
import numpy as np
import gc

file_length = 200000
vec_dim = 2048

FILE_PATH = '/data/workspace/test'
OUT = '/data/workspace/out_test'
OUT_SMILES = 'out_smiles'
OUT_IDS = 'out_ids'
OUT_NPY = 'out_npy'


def thread_runner(smiles, ids, filename):
    thread_num = math.ceil(len(smiles)/file_length)
    with ProcessPoolExecutor(thread_num) as executor:
        for i in range(thread_num):
            try:
                executor.submit(get_smiles_fp, smiles[i*file_length:i*file_length+file_length], ids[i*file_length:i*file_length+file_length], filename, i)
            except:
                executor.submit(get_smiles_fp, smiles[i*file_length:len(smiles)], ids[i*file_length:len(ids)], filename, i)


def get_smiles_fp(smiles, ids, filename, num):
    hex_fps = []
    new_ids = []
    new_smiles = []
    count = 0
    for smile in smiles:
        m = Chem.MolFromSmiles(smile)
        if m is None:
            count = count + 1
            continue
        # the method for generating fingerprints - morgan_fp or rdk_fp
        # fp2 = AllChem.GetMorganFingerprintAsBitVect(m, 2, vec_dim)
        fp2 = Chem.RDKFingerprint(m, fpSize=vec_dim)
        hex_fp = DataStructs.BitVectToFPSText(fp2)
        hex_fps.append(hex_fp)
        new_ids.append(ids[count])
        new_smiles.append(smile)
        count = count + 1
    hex_fps = np.array(hex_fps)
    np.save(OUT + '/' + OUT_NPY + '/' + filename + "%03d"%num +'.npy', hex_fps)
    save_file(new_smiles, OUT + '/' + OUT_SMILES + '/' + filename + "%03d"%num +'.smi')
    save_file(new_ids, OUT + '/' + OUT_IDS + '/' + filename + "%03d"%num +'.txt')
    del hex_fps
    del new_smiles
    del new_ids
    gc.collect()


def save_file(news, fname):
    with open(fname,'a') as f:
        for new in news:
            f.write(new + '\n')


def get_smiles(path, fname):
    smiles = []
    ids = []
    filename = path + '/' + fname
    with open(filename, "r") as infile:
        for line in infile:
            parts = line.split()
            try:
                ids.append(parts[1])
                smiles.append(parts[0])
            except :
                print(filename)
    return smiles, ids


def get_files_fp(file_path):
    filenames = os.listdir(file_path)
    filenames.sort()
    for filename in filenames:
        smiles, ids = get_smiles(file_path,filename)
        # print(len(smiles), len(ids), filename[0])
        thread_runner(smiles, ids, filename[0])


def main():
    if not os.path.exists(OUT):
        os.mkdir(OUT)
        os.mkdir(OUT + '/' + OUT_NPY)
        os.mkdir(OUT + '/' + OUT_SMILES)
        os.mkdir(OUT + '/' + OUT_IDS)
        
    time1 = time.time()
    get_files_fp(FILE_PATH)
    time2 = time.time()
    print("total time: ", time2 - time1)


if __name__ == "__main__":
    main()
