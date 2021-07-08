from random import shuffle
from itertools import permutations, product
from pandas import DataFrame
import numpy as np

#Skip function 
def skipper(pssm, skip):
    skips_fragments = []
    skip_steps = list(permutations(pssm.index, skip))
    skip_steps = list(map(list, skip_steps))
    [com.sort() for com in skip_steps]
    skip_steps = list(map(list, list(set(map(tuple, skip_steps)))))

    for step in skip_steps:

        pssm_step = pssm
        if skip != 0:
            for skip_module in step:

                pssm_step = pssm_step[pssm_step.index != skip_module]
        skips_fragments.append(pssm_step)
        
    return skips_fragments

#Making combinations of modules
def create_variants(original_seq, len_place,):
    
    seq = []
    
    for s in original_seq:
        if len(s) != 0:
            seq.append(s)
            
    N = sum([len(s) for s in seq])
    N_cnt = len(seq) + 1
    orders = list(
        permutations([i for i in range(len(seq))])
    )
    free_slots = len_place - N
    
    variants = []
    for rep in product(
        [i for i in range(free_slots + 1)],
        repeat = N_cnt
    ):
        if sum(rep) == free_slots:
            variants.append(rep)
    concatenates = []
    
    for v in variants:
        for o in orders:
            concat = []
            for j in range(len(v)-1):
                concat += ['nan'] * v[j]
                concat += seq[o[j]]
            concat += ['nan'] * v[-1]
            concatenates.append(concat)
    
    
    
    concatenates = list(set (map(tuple, concatenates)))
    concatenates.sort()
    return concatenates

    variants = set(variants)

    return variants

#Standartization of monomers
def make_standard(sequences, subtrate_stack):
    
    check_1 = 0
    check_2 = 0
    for seq in sequences:
        for sub_AS in subtrate_stack:
            for sub in seq:

                if sub_AS in sub:
                    
                    sequences[sequences.index(seq)][seq.index(sub)] =  sub_AS

    return sequences

def make_combine(sequences, subtrate_stack, pssm='None', delta=3):
    for var in sequences:

        sequences[var] = make_standard(sequences[var], subtrate_stack)
    
    N_refers = []

    for var in sequences:
        
        N_refer = 0

        for cont in sequences[var]:
            for mon in cont:
                
                N_refer += 1
                
                if mon not in subtrate_stack:

                    cont[cont.index(mon)] = 'nan'

        N_refers.append(N_refer)
    
    if type(pssm) == DataFrame:
        if min(N_refers) + delta < len(pssm):
            return None

        else:

            answ = []

            for var in sequences:
                
                nx = []
                
                for i in sequences[var]:
                    
                    nx.append(i)

                answ.extend(create_variants(nx, len(pssm)))
            answ = list(set(answ))
            return answ
        
    else:
        
        return  min(N_refers)
    
#Making shuffle funtions
def shuffle_matrix(pssm_profile):
    
    profile = pssm_profile.copy()

    while np.min(profile.values == pssm_profile.values) == True:
        for idx in profile.index:

            row_vals = list(profile.iloc[idx].values[1: ])
            shuffle(row_vals)
            row = {k : v for k, v in zip(profile.iloc[idx].keys()[1: ], row_vals)}
            
            for sub in profile.keys()[1: ]:

                profile[sub][idx] = row[sub]

    return profile
        
def pssm(seq, pssm_profile):

    TS = []
    cut = 0
    N = pssm_profile.shape[0]
    #Calculating C-score(cluster-score)
    
    while N - cut - len(seq) != -1:
        seq_cnt = 0
        target_sum = 0

        for line in pssm_profile.index:
            if line < cut:
                continue
            try:
                
                target_sum += pssm_profile[seq[seq_cnt]][line]
                
            except KeyError:
                
                
                target_sum += 0
                
            except IndexError:
                
                target_sum += 0

            seq_cnt += 1
        
        cut += 1
        TS.append(target_sum)
    
    return max(TS)
    