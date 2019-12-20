import sys
import os
import pickle
import math
import multiprocessing as mp
from rdkit import Chem
from rdkit.Chem import AllChem

def score(chunk_path, db_shingles, rooted, radius, weighted, cut_off, sender):
    res_path = '_res'+chunk_path
    radius_constr = radius + 1

    with open(chunk_path, 'r') as chunk, open(res_path, 'w') as results:
        for line in chunk:

            qry_shingles = set()
            smi = line.split('\t')[0].rstrip()
            mol = Chem.MolFromSmiles(smi)
            if mol:
                for atm_idx in range(mol.GetNumAtoms()):
                    for N in range(1,radius_constr):
                        bonds = AllChem.FindAtomEnvironmentOfRadiusN(mol, N, atm_idx)

                        if not bonds:
                            break

                        # the reportedly faster method
                        atoms = set()
                        for bond_id in bonds:
                            bond = mol.GetBondWithIdx(bond_id)
                            atoms.add(bond.GetBeginAtomIdx())
                            atoms.add(bond.GetEndAtomIdx())
                        
                        if rooted:
                            new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol,list(atoms),bonds,0,0,False,False,atm_idx,True,False,False)
                        else:
                            new_shingle = Chem.rdmolfiles.MolFragmentToSmiles(mol,list(atoms),bonds,0,0,False,False,-1,True,False,False)
                        
                        '''
                        # the slower method
                        atm_map = {}
                        submol = Chem.PathToSubmol(mol,bonds,atomMap=atm_map)
                        
                        new_shingle = Chem.MolToSmiles(submol,rootedAtAtom=atm_map[atm_idx],canonical=True)
                        '''

                        qry_shingles.add(new_shingle)
                    
                # calculate shingle count averaged score            
                avg_score = 0
                if qry_shingles:
                    sum_scores = 0
                    # using log10 of shingle frequency
                    if weighted:
                        for shingle in qry_shingles:
                            # if key not present, add 0 per default
                            sum_scores += db_shingles.get(shingle, 0)
                    # working binary (i.e. if present -> count++ )
                    else:
                        for shingle in qry_shingles:
                            if shingle in db_shingles:
                                sum_scores += 1
                    avg_score = sum_scores/len(qry_shingles)

                if cut_off == 0.0 or cut_off <= avg_score:
                    results.write(line.rstrip()+'\t'+str(avg_score)+'\n')
            
    sender.send(res_path)


if __name__ == "__main__":
    p_in = ''
    p_out = ''
    radius = 3          # default
    rooted = True       # default
    weighted = True     # default
    min_freq = True     # default
    cut_off = 0.0       # default
    
    # parse arguments
    arglen = len(sys.argv)
    
    for i, arg in enumerate(sys.argv[1:],1):
        last_inicates_explicit = False
            
        if arg.startswith('--'):
            if arg.startswith('--radius='):
                radius = int(arg.replace('--radius=',''))
                if radius <= 1:
                    print('Please set radius >= 2 and <= shingle database maximum radius...')
                    exit()
            elif arg.startswith('--rooted=') and arg.endswith('alse'):
                rooted = False
            elif arg.startswith('--weighted=') and arg.endswith('alse'):
                weighted = False
            elif arg.startswith('--considerRareShingles=') and arg.endswith('rue'):
                min_freq = False
            elif arg.startswith('--cutOffScore='):
                cut_off = float(arg.split('--cutOffScore=')[-1])

        elif arg.startswith('-'):
            last_inicates_explicit = True

            if arg == '-i' and  i < arglen and os.path.isfile(sys.argv[i+1]):
                p_in = sys.argv[i+1]

            elif arg == '-o' and  i < arglen and sys.argv[i+1] != p_in:
                p_out = sys.argv[i+1]

        elif not last_inicates_explicit and p_in == '' and os.path.isfile(sys.argv[i]):
            p_in = sys.argv[i]
            last_inicates_explicit = False

        elif not last_inicates_explicit and p_out == '' and sys.argv[i] != p_in:
            p_out = sys.argv[i]
            last_inicates_explicit = False
  
    if p_in == '':
        print('Please give source file path as argument...')
        exit()
    if p_out == '':
        p_out = '.'.join(p_in.split('.')[:-1])+'_res.'+p_in.split('.')[-1]
        print('No path given to write scores. Will write to:\n'+p_out)

    # Load respective shingle db
    db_shingles = False
    if rooted:
        if min_freq:
            with open('shingle_libs/chembl_24_1_shingle_scores_log10_rooted_nchir_min_freq_100.pkl', 'rb') as pyc:
                db_shingles = pickle.load(pyc)
        else:      
            #with open('UNPD_shingleDB_rt_r3.pyc', 'rb') as pyc:
            with open('shingle_libs/chembl_24_1_shingle_scores_log10_rooted_nchir.pkl', 'rb') as pyc:
                db_shingles = pickle.load(pyc)
    else:
        with open('shingle_libs/chembl_24_1_shingle_scores_log10_nrooted_nchir.pkl', 'rb') as pyc:
            db_shingles = pickle.load(pyc)

    if not db_shingles:
        print('Could not read pickled shingle db...')
        exit()

    print('Database involves '+str(len(db_shingles))+' shingles.\nRadius will be: '+str(radius)+'\nUsing rooted: '+str(rooted)+'\nUsing weighted score: '+str(weighted))
    
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
        p = mp.Process(target=score, args=(path, db_shingles, rooted, radius, weighted, cut_off, sender))
        jobs.append(p)
        pipes.append(reciever)
        p.start()
    for proc in jobs:
        proc.join()
    results = [c.recv() for c in pipes]


    # merge results to new file and remove temp result chunks
    p_in_split = p_in.split('.')

    with open(p_out, 'w') as orig:
        for path in results: 
            with open(path, 'r') as fi_res:
                for line in fi_res:
                    orig.write(line)
            os.remove(path)
        results= []

    # remove temp source files
    for path in temp_paths:
        os.remove(path)
