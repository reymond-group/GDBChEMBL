import sys
import operator
from rdkit import Chem
from rdkit.Chem import AllChem
import math
import pickle

p_in = sys.argv[1]
p_out = sys.argv[2]

db_shingles = {}
sh_count = 0
with open(p_in, 'r') as fi_in:
    #fi_in.readline() # header
    for i, line in enumerate(fi_in):
        smi = line.split('\t')[0].rstrip()
        mol = Chem.MolFromSmiles(smi)

        if mol:
            for atm in Chem.CanonicalRankAtoms(mol):
                for N in range(1,4):
                    bonds = AllChem.FindAtomEnvironmentOfRadiusN(mol, N, atm)

                    if not bonds:
                        break
                    
                    # the faster method... 
                    atoms = set()
                    for bond_id in bonds:
                        bond = mol.GetBondWithIdx(bond_id)
                        atoms.add(bond.GetBeginAtomIdx())
                        atoms.add(bond.GetEndAtomIdx())
                    shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol,list(atoms),bonds,0,0,False,False,atm,True,False,False)
                    '''
                    # the slower method
                    submol = Chem.PathToSubmol(mol,bonds)
                    shingle = Chem.MolToSmiles(submol, False,False,atm,True,False,False)
                    '''
                    sh_count += 1

                    if shingle in db_shingles.keys():
                        db_shingles[shingle] += 1
                    else:
                        db_shingles[shingle] = 1

for key, value in db_shingles.items():
    db_shingles[key] = math.log(value, 10)+1 # default base: e

with open(p_out, 'wb') as pyc:
    pickle.dump(db_shingles,pyc)
