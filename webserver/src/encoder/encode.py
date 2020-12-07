import os
import numpy as np
from diskcache import Cache
from rdkit import DataStructs
from rdkit import Chem
import math
from rdkit.Chem import AllChem
from common.config import VECTOR_DIMENSION


def smiles_to_vec(smiles):
    mols = Chem.MolFromSmiles(smiles)
    # fp = AllChem.GetMorganFingerprintAsBitVect(mols, 2, VECTOR_DIMENSION)
    fp = Chem.RDKFingerprint(mols, fpSize=VECTOR_DIMENSION)
    hex_fp = DataStructs.BitVectToFPSText(fp)
    # print(hex_fp)
    vec = bytes.fromhex(hex_fp)
    return vec