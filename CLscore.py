import sys
import os
import pickle
import math
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import AllChem

def score(chunk_path, db_shingles, sender):
    res_path = '_res'+chunk_path
    with open(chunk_path, 'r') as chunk, open(res_path, 'w') as results:
        for line in chunk:
            qry_shingles = []
            smi = line.split('\t')[0].rstrip()
            mol = Chem.MolFromSmiles(smi)
            for atm in Chem.CanonicalRankAtoms(mol):
                for N in range(1,4):
                    bonds = AllChem.FindAtomEnvironmentOfRadiusN(mol, N, atm)

                    if not bonds:
                        break

                    atoms = set()
                    for bond_id in bonds:
                        bond = mol.GetBondWithIdx(bond_id)
                        atoms.add(bond.GetBeginAtomIdx())
                        atoms.add(bond.GetEndAtomIdx())
                    new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol,list(atoms),bonds,0,0,False,False,atm,True,False,False)
                    '''
                    submol = Chem.PathToSubmol(mol,bonds)
                    new_shingle = Chem.MolToSmiles(submol) # UNROOTED!!!!
                    '''
                    if new_shingle not in qry_shingles:
                        qry_shingles.append(new_shingle)
            
            avg_score = 0
            if qry_shingles:
                sum_scores = 0
                for shingle in qry_shingles:
                    # if key not present, add 0 per default
                    sum_scores += db_shingles.get(shingle, 0)
                avg_score = sum_scores/len(qry_shingles)

            results.write(line.rstrip()+'\t'+str(avg_score)+'\n')
            
    sender.send(res_path)


if __name__ == "__main__":
    db_shingles = False
    with open('chembl_24_1_shingle_scores_log10_rooted_nchir.pyc', 'rb') as pyc:
        db_shingles = pickle.load(pyc)

    print(str(len(db_shingles)))
    if not db_shingles:
        print('Could not read pickled shingle db...')
        exit()

    p_in = sys.argv[1]

    length_in = 0
    with open(p_in, 'r') as fi_in: 
        for line in fi_in:
            length_in += 1

    cpu_count = mp.cpu_count()
    chunk_len = math.ceil(length_in / cpu_count)

    temp_paths = []
    with open(p_in, 'r') as fi_in:  
        #fi_in.readline() # header TEMP  
        for i in range(1,cpu_count):
            p_temp = '_temp'+str(i)+'_'+p_in.split('/')[-1]
            temp_paths.append(p_temp)
            with open(p_temp, 'w') as temp:
                for j, line in enumerate(fi_in):
                    temp.write(line)
                    if j >= chunk_len:
                        break

    # calculate scores in parallel
    jobs = []
    pipes = []
    for path in temp_paths:
        reciever, sender = mp.Pipe(False)
        p = mp.Process(target=score, args=(path, db_shingles, sender))
        jobs.append(p)
        pipes.append(reciever)
        p.start()
    for proc in jobs:
        proc.join()
    results = [c.recv() for c in pipes]


    # merge results to new file and remove temp result chunks
    p_in_split = p_in.split('.')
    p_res = '.'.join(p_in_split[:-1])+'_res.'+p_in.split('.')[-1]
    with open(p_res, 'w') as orig:
        for path in results: 
            with open(path, 'r') as fi_res:
                for line in fi_res:
                    orig.write(line)
            os.remove(path)
        results= []

    # remove temp source files
    for path in temp_paths:
        os.remove(path)