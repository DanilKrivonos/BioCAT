import numpy as np
import logging
from subprocess import call 
import os
from os import listdir
from pandas import read_csv
from Bio.SearchIO import parse

a_logger = logging.getLogger()
a_logger.setLevel(logging.DEBUG)

""" The main role of this functions is generate PSSM """

# The function deletes protocluster, which is duplicate of candidate cluster
def check_cluster_duplicates(pssm_path):
    """
    The function remove duplicates of PSSMs.

    Parameters
    ----------
    pssm_path : str
        Path to PSSMs.
    """
    pssms = listdir(pssm_path)
    
    for pssm in pssms:
        if 'cand' not in pssm:
            continue

        CC_pssm = read_csv(pssm_path + pssm, sep='\t')
        
        if len(pssm.split('_')) == 5:
        
            protoID = pssm.split('_')[4].split('.')[0]
            
        else:
            
            protoID = pssm.split('_')[4]

        for p_pssm in pssms:
            if 'proto_{}_'.format(protoID) not in p_pssm:
                continue
            try:
                
                proto_pssm = read_csv(pssm_path + p_pssm, sep='\t')
            
            except FileNotFoundError: # if file was removed
                continue
            try:
                if all(CC_pssm == proto_pssm) == True:
                    
                    a_logger.debug(p_pssm + ' is duplicate')
                    call('rm ' + pssm_path + p_pssm, shell=True)

            except ValueError: # if the cluster was splited, its will differ
                continue
        

def PSSM_make(search, aminochain, out, delta, substance_name):
    """
    The function genetates PSSM profiles from hmm otputs.

    Parameters
    ----------
    search : str
        Path to HMM profile directory
    aminochain : str
        Length of minimal possible core peptide chain variant for the substances. 
    out : str
        Path to BioCAT output directory.
    delta : int
        Differense between length of core peptide chain and size of cluster.
    substance_name : str
        Name of analyzing substance.
    """
    hmms = os.path.dirname(os.path.abspath(__file__)) + '/../HMM/Bacteria_HMM'
    neg = os.path.dirname(os.path.abspath(__file__)) + '/../HMM/Bac_negative.tsv'
    negative = read_csv(neg, sep='\t')
    output = out
    path = out + '/table.tsv'
    substr_hmm = listdir('{}'.format(hmms))
    substr = []
    # Making substrates list
    for subst in substr_hmm:

        substr.append(subst.split('_')[1][: -4])

    table = read_csv(path, sep='\t')
    subst_eval = {}
    
    # Making dict with all BGCs
    for id_new in table[table['Domain name'].str.contains('AMP-binding')].index:
        if table['ID'][id_new] not in subst_eval:

            subst_eval[table['ID'][id_new]] = {} 
        
        key = table['Domain name'][id_new]
        subst_eval[table['ID'][id_new]][key] = {}

        for sub in substr:

            subst_eval[table['ID'][id_new]][key][sub] = []
    #This dict contains only possible BGCs, which length bigger or equel amino chain of molecule
    cop =  subst_eval.copy()
    
    for s in  cop:
        if len(cop[s]) < aminochain:

            subst_eval.pop(s)
            continue
        
        if len(cop[s]) > aminochain + delta:
            
            subst_eval.pop(s)
            continue       
    
    # Generator of ratio
    a_logger.debug('Calculation ratio scores for PSSM ...')

    for sub in substr:
        
        TN = [-float(i) for i in list(negative[negative.Substrate == sub].List_of_negative)[0].split(',')]
        TN = np.array(TN)

        for substance in subst_eval.keys():
            for module in  subst_eval[substance].keys():

                open_hmm  = parse('{}/{}.out'.format(search, sub), 'hmmer3-text')
                
                for i in open_hmm:
                    for t in i:
                        if module == t.hsps[0].hit.id:
                            
                            TP = -np.log(t.evalue)
                
                ratio_metric = len(TN[TN <= TP])/ len(TN)
                subst_eval[substance][module][sub] = ratio_metric

    a_logger.debug('Recording PSSM ...\n')
    
    #Making directory for PSSMs
    PSSMs_out = output + '/PSSM_{}/'.format(substance_name)
    try:

        os.mkdir('{}'.format(PSSMs_out))

    except FileExistsError:

        a_logger.debug('The output directory already exists')
    #Recording PSSMs
    for substance in subst_eval.keys():
        with open('{}/PSSM_A_{}.csv'.format(PSSMs_out, substance), 'w') as test:

            test.write('Module name\t')
            test.write('{}\n'.format('\t'.join(substr)))

            for module in  subst_eval[substance].keys():
                
                test.write('{}'.format(module))
                
                for sub in substr: 

                    test.write('\t')
                    test.write('{}'.format(subst_eval[substance][module][sub]))
                    
                test.write('\n')
    a_logger.debug('Checking cluster duplicates ...')
    check_cluster_duplicates(PSSMs_out)# delete duplicate 
    a_logger.debug('All PSSMs were recorded!')
