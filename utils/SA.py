import os
import glob
import gzip
import math
import pickle
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import yaml

_fscores = None

def readFragmentScores(name='fpscores'):
    global _fscores
    if name == "fpscores":
        name = os.path.join(os.getcwd(), name)
    data = pickle.load(gzip.open('%s.pkl.gz' % name, 'rb'))
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict

def numBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro

def calculateSA(smiles):
    if _fscores is None:
        readFragmentScores()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0

    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += _fscores.get(sfp, -4) * v
    score1 /= nf

    nAtoms = mol.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    ri = mol.GetRingInfo()
    nBridgeheads, nSpiro = numBridgeheadsAndSpiro(mol, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = 0.

    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    score3 = 0.
    if nAtoms > len(fps):
        score3 = math.log(float(nAtoms) / len(fps)) * .5

    sascore = score1 + score2 + score3

    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.

    if sascore > 8.:
        sascore = 8. + math.log(sascore + 1. - 9.)
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    return sascore

def load_config(config_path):
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def calculate_SA(config):
    config_output_path = config['output']['output_path']
    input_path = os.path.join(config_output_path, 'General-molecular-quality-metrics')
    uniqueness_base_folder = os.path.join(input_path, 'Uniqueness')
    output_base_folder = os.path.join(input_path, 'SA')
    os.makedirs(output_base_folder, exist_ok=True)

    for i in range(32):
        numeric_subfolder = str(i)
        uniqueness_folder = os.path.join(uniqueness_base_folder, numeric_subfolder)
        output_folder = os.path.join(output_base_folder, numeric_subfolder)

        if not os.path.exists(uniqueness_folder):
            print(f" {uniqueness_folder} non-existent. skip")
            continue
        
        os.makedirs(output_folder, exist_ok=True)

        all_sa_scores_file_path = os.path.join(output_folder, 'SA.txt')
        
        with open(all_sa_scores_file_path, 'w') as f_all_sa_scores:
            txt_files = glob.glob(os.path.join(uniqueness_folder, '*.txt'))
            
            if not txt_files:
                print(f" {uniqueness_folder} without .txt . skip")
                continue

            for txt_file_path in txt_files:
                with open(txt_file_path, 'r') as txt_file:
                    lines = txt_file.readlines()
                    
                    for line in lines:
                        data = line.strip().split('\t')
                        if len(data) >= 2:
                            identifier, smiles = data[0], data[1]  
                            sa_score = calculateSA(smiles)
                            if sa_score is not None:
                                f_all_sa_scores.write(f"{identifier}\tSMILES:\t{smiles}\tSA Score:\t{sa_score}\n")
                            else:
                                print(f"Error processing SMILES: {smiles}")
                        else:
                            print(f"Invalid data format in file: {txt_file_path}")
                
        #print(f" {numeric_subfolder} Processing completed。\n")


