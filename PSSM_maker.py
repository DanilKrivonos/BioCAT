import numpy as np
from subprocess import call 
import os
from os import listdir
from pandas import read_csv
from Bio.SearchIO import parse

def PSSM_make(search, aminochain, out, delta, hmms='./HMM/', neg='./negative.tsv'):

    negative = read_csv(neg, sep='\t')
    output = out
    path = out + 'table.tsv'
    substr = listdir('{}'.format(hmms))
    table = read_csv(path, sep='\t')
    subst_eval = {}

    #making dict with all BGCs
    for id_new in table[table['Domain name'].str.contains('AMP-binding')].index:
        
        if table['ID'][id_new] not in subst_eval:
            subst_eval[table['ID'][id_new]] = {} 
            
        key = '{}_{}'.format(table['Domain name'][id_new], table['Single aa prediction'][id_new])

        if key in subst_eval:
            continue

        subst_eval[table['ID'][id_new]][key] = {}
        
        for sub in substr:
            subst_eval[table['ID'][id_new]][key][sub] = []
    #This dict only with possible BGCs, which length bigger or equel amino chain of molecule
    subst_eval_new = subst_eval.copy()
    for s in  subst_eval_new:
        if len(subst_eval_new[s]) < aminochain:

            subst_eval.pop(s)
            continue
        
        if len(subst_eval[s]) > aminochain + delta:

            subst_eval.pop(s)
            continue
    
    #Generator of n-score
    print('Calculation n-scores ...')
    for sub in substr:
        
        FP = [-float(i) for i in negative.loc[negative['Substrate'] == '{}'.format(sub)]['List_of_negative'].values[0].split('_')]
        
        for substance in subst_eval.keys():
            print('Calculation n-scores for {}\n'.format(substance))
            for module in  subst_eval[substance].keys():
                
                open_hmm  = parse('{}/{}vs{}.out'.format(search, sub, module, sub), 'hmmer3-text')

                for i in open_hmm:
                    for t in i:

                        TP = -np.log(t.evalue)

                FP = np.array(FP)
                ratio_metric = len(FP[FP <= TP])/ len(FP)
                subst_eval[substance][module][sub] = ratio_metric

    print('Recording PSSM ...\n')

    #Making directory for PSSMs
    PSSMs_out = output + 'PSSM/'
    try:

        os.mkdir('{}'.format(PSSMs_out))

    except FileExistsError:

        print('The output directory already exists')
    #Recording PSSMs
    for substance in subst_eval.keys():
        print('For {}'.format(substance))
        with open('{}/PSSM_A_{}.csv'.format(PSSMs_out, substance), 'w') as test:

            test.write('Module name\t')
            test.write('{}\n'.format('\t'.join(substr)))

            for module in  subst_eval[substance].keys():

                test.write('{}'.format(module))
                
                for sub in substr: 

                    test.write('\t')
                    test.write('{}'.format(subst_eval[substance][module][sub]))
                    
                test.write('\n')
                                
    print('All PSSM were recorded sucessfully!')