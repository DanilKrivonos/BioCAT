from random import shuffle
from itertools import permutations, product
from pandas import DataFrame
from os import listdir
from numpy import linspace
import numpy as np
import os
from multiprocessing import Process, Manager

#Skip mode function 
def skipper(pssm, skip):
    """
    Skip function. In the rare cases, then one or more of modules in cluster 
    is skipped (https://doi.org/10.1186/s12934-018-0929-4).

    Parameters
    ----------
    pssm : pandas DataFrame
        PSSM profile.
    skip : int
        Number of presumptive skip.
    Returns
    -------
    skips_fragments : list
        Each possible variants of pssm with presumptive number of skip.
    """
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
def create_variants(original_seq, len_place):
    """
    Combined fragments to possible sequences.

    Parameters
    ----------
    original_seq : list
        Sequence of possible fragments.
    len_place : int
        Length of PSSM profile (length of cluster).
    Returns
    -------
    concatenates : list
        Possible combined sequence from continual fragments.
    """
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

#Making possible combinations for potential clusters
def make_combine(sequence, length_max, pssm, delta=3):
    """
    Create possible combinations of core peptide sequences.

    Parameters
    ----------
    sequence : list
        Sequences of one of biosynthesis type (e.g. A, B, or C).
    length_max : int
        Length of PSSM profile (length of cluster).
    pssm : pandas DataFrame
        PSSM profile.
    delta : int
        Differense between length of core peptide chain and size of cluster.
    Returns
    -------
    answ : list
        Possible combined sequence from continual fragments.
    """
    if length_max + delta < len(pssm): #if minimal possible length of peptide sequence much less, then cluster 
        return None

    else:

        answ = []
        
        for var in sequence:
            
            nx = []
            
            for i in sequence[var]:
                
                nx.append(i)
            
            answ.extend(create_variants(nx, len(pssm)))

        answ = list(set(answ))

    return answ #Returen combinations

#Getting maximum size of putative sequence from variant of biosynthesis type
def get_max_aminochain(bios_path):
    """
    The function return length of minimum fragmnt of core peptide chain.

    Parameters
    ----------
    bios_path : list
        Sequences of one of biosynthesis type (e.g. A, B, or C).
    Returns
    -------
    max(lens_of_varints) : int
        Maximum possible lenght from possible variants.
    """
    lens_of_varints = []

    for var in bios_path:
        
        len_of_varint_seq = 0

        for cont in bios_path[var]:

            len_of_varint_seq += len(cont)
        
        lens_of_varints.append(len_of_varint_seq)

    return  max(lens_of_varints) # Return maximum lenght from possible variants
#Getting minimla size of putative sequence
def get_minim_aminochain(PeptideSeq):
    """
    The function return length of minimum fragment of core peptide chain.

    Parameters
    ----------
    PeptideSeq : dict
        Core peptide chains for different biosynthesis types (e.g. A, B, or C).
    Returns
    -------
    min(lens_of_varints) : int
        Minimal possible lenght from possible variants.
    """
    lens_of_varints = []

    for bios_path in PeptideSeq:
        for var in PeptideSeq[bios_path]:
            
            len_of_varint_seq = 0

            for cont in PeptideSeq[bios_path][var]:

                len_of_varint_seq += len(cont)

            lens_of_varints.append(len_of_varint_seq)

    return  min(lens_of_varints) # Return minimal lenght from possible variants
#SHUFFLUNG FUNCTION
def shuffle_matrix(pssm_profile, ShufflingType):
    """
    The function return shuffled PSSM profile.

    Parameters
    ----------
    pssm_profile : pandas DataFrame
        PSSM profile.
    ShufflingType : str
        Variant of shuffling.
    Returns
    -------
    profile : pandas DataFrame
        Shuffled PSSM profile.
    """
    # Module shuffling
    if ShufflingType == 'module':
        
        profile = pssm_profile.copy()
        cols = profile.columns[1: ]
        
        while np.min(profile.values == pssm_profile.values) == True:
            for col in cols:

                col_vals = list(profile[col])
                shuffle(col_vals)
                profile[col] = col_vals
            if len(profile) < 3:
                break
    # Substrate shuffling
    elif ShufflingType == 'substrate':

        profile = pssm_profile.copy()

        while np.min(profile.values == pssm_profile.values) == True:
            for idx in pssm_profile.index:

                row_vals = list(profile.iloc[idx].values[1: ])
                shuffle(row_vals)
                row = {k : v for k, v in zip(profile.iloc[idx].keys()[1: ], row_vals)}

                for sub in profile.keys()[1: ]:

                    profile.loc[idx, sub] = row[sub]

    return profile 
#
def get_score(seq, pssm_profile, type_value):
    """
    Alighning peptide sequence on PSSM profile.

    Parameters
    ----------
    seq : str
        Core peptide chain.
    pssm_profile : str
        PSSM profile.
    type_value : str or None
        Variant of score modification (log or default).
    Returns
    -------
    profile : pandas DataFrame
        Shuffled PSSM profile.
    """
    target_sum = 0
    seq_cnt = 0

    for line in pssm_profile.index:
        
        try:
            if type_value == 'log':
                if pssm_profile[seq[seq_cnt]][line] != 0:

                    target_sum += np.log(pssm_profile[seq[seq_cnt]][line])

                else:

                    target_sum += -np.inf
                
            elif type_value == None:

                target_sum += pssm_profile[seq[seq_cnt]][line]
            
        except:
            
            target_sum += 0

        seq_cnt += 1
      
    return target_sum

# Multple treading generating shuflling mtrixes
def get_shuffled_matrix(pssm_mat, iterations, return_dict, ShufflingType):
    """
    The functuion generate massive of shuffled matrix.

    Parameters
    ----------
    pssm_mat : pandas DataFrame
        PSSM profile.
    iterations : int
        Number of iterations of shuffling.
    return_dict : dict
        Manager of multitreating.
    ShufflingType : str
        Variant of shuffling.
    """
    shuffled_matrix = []
    for i in range(iterations):

        shuffled_matrix.append(shuffle_matrix(pssm_mat, ShufflingType))

    return_dict[os.getpid()] = shuffled_matrix

def multi_thread_shuffling(pssm_mat, ShufflingType, iterations=100, threads=1):
    """
    The functuion parallel matrix shuffling.

    Parameters
    ----------
    pssm_mat : pandas DataFrame
        PSSM profile.
    ShufflingType : str
        Variant of shuffling.
    iterations : int
        Number of iterations of shuffling.
    threads : int
        Number of threads used.
    Returns
    -------
    shuffled_matrix : list
        List of shuffled matrix.
    """
    procs = []
    manager = Manager()
    return_dict = manager.dict()
    thread_iterations = [int(iterations/threads)] * threads    

    thread_iterations[0] += iterations - sum(thread_iterations)

    for i in thread_iterations:
        proc = Process(target=get_shuffled_matrix, args=(pssm_mat, i, return_dict, ShufflingType))
        procs.append(proc)
        proc.start()
        
    for proc in procs:
        proc.join()
        
    shuffled_matrix = []

    for v in return_dict.values():
        shuffled_matrix += v
    
    return shuffled_matrix

# Multple treading calculation scores for shuffled matrix
def get_shuffled_scores(MaxSeq, return_dict, ShuffledMatrix, type_value):
    """
    The functuion calculate scores for shuffled matrix.

    Parameters
    ----------
    MaxSeq : list
        Suquence with maximum possible value.
    return_dict : dict
        Manager of multitreating.
    ShuffledMatrix  : pandas DataFrame
        Variant of shuffled matrix.
    type_value : str or None
        Variant of score modification (log or default).
    """
    shuffled_scores = []
    for i in range(len(ShuffledMatrix)):
        shuffled_scores.append(get_score(MaxSeq, ShuffledMatrix[i], type_value))

    shuffled_scores = np.array(shuffled_scores)
    return_dict[os.getpid()] = shuffled_scores


def multi_thread_calculating_scores(MaxSeq, ShuffledMatrix, type_value, iterations=100, threads=1):
    """
    The functuion parallel calculation of scores for shuffled matrix.

    Parameters
    ----------
    MaxSeq : list
        Suquence with maximum possible value.
    ShuffledMatrix : pandas DataFrame
        Variant of shuffled matrix.
    type_value : str or None
        Variant of score modification (log or default).
    iterations : int
        Number of iterations of shuffling.
    threads : int
        Number of threads used.
    Returns
    -------
    shuffled_scores : list
        List of scores for shuffled matrix.
    """
    procs = []
    manager = Manager()
    return_dict = manager.dict()
    #thread_iterations = [int((len(ShuffledMatrix))/threads)] * threads    
    batches = linspace(0, len(ShuffledMatrix), threads)
    #thread_iterations[0] += iterations - sum(thread_iterations)

    edges = [
        (batches[i], batches[i+1]) for i in range(len(batches) - 1)
    ]
    for i in range(len(edges)):

        proc = Process(target=get_shuffled_scores, args=(MaxSeq, return_dict, ShuffledMatrix[int(edges[i][0]): int(edges[i][1])], type_value))
        procs.append(proc)
        proc.start()
        
    for proc in procs:
        proc.join()

    #print(len(return_dict))
    shuffled_scores = np.concatenate([return_dict[k] for k in return_dict.keys()])
    
    return shuffled_scores