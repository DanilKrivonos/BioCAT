import numpy as np
from subprocess import call 
import os
from os import listdir
from pandas import read_csv
from Bio.SearchIO import parse

def PSSM_make(search, aminochain, out, delta, hmms='./HMM/', neg='./negative.tsv'):

    negative = read_csv(neg, sep='\t')
    output = out
    path = out + '/table.tsv'
    substr = listdir('{}'.format(hmms))
    table = read_csv(path, sep='\t')
    subst_eval = {}
    
    #making dict with all BGCs
    for id_new in table[table['Domain name'].str.contains('AMP-binding')].index:
        if table['ID'][id_new] not in subst_eval:

            subst_eval[table['ID'][id_new]] = {} 
        
        key = '{}_{}'.format(table['Domain name'][id_new], table['Single aa prediction'][id_new])
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
    
    #Generator of n-score
    print('Calculation scores for PSSM ...')

    for sub in substr:
        
        FP = [-float(i) for i in list(negative[negative.Substrate == sub].List_of_negative)[0].split(',')]

        for substance in subst_eval.keys():
            for module in  subst_eval[substance].keys():
                
                open_hmm  = parse('{}/{}.out'.format(search, sub), 'hmmer3-text')
            
                for i in open_hmm:
                    for t in i:
                        if t.hsps[0].hit.id == '{}_{}'.format(substance, module):
                       
                            TP = -np.log(t.evalue)
                            break
                
                FP = np.array(FP)
                ratio_metric = len(FP[FP <= TP])/ len(FP)
                subst_eval[substance][module][sub] = ratio_metric

    print('Recording PSSM ...\n')
    
    #Making directory for PSSMs
    PSSMs_out = output + '/PSSM/'
    try:

        os.mkdir('{}'.format(PSSMs_out))

    except FileExistsError:

        print('The output directory is already exists')
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
