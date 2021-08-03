import argparse
from pandas import read_csv
from subprocess import call
from os import listdir, mkdir

def HMM_make(path, output, cpu, hmms='./HMM/'):

    table = read_csv(path + '/table.tsv', sep='\t')
    #making work directory to fasta files
    substrates = []
        
    hmm_out = '{}/HMM_results/'.format(output)
    try:
        
        mkdir('{}'.format(hmm_out))
        
    except:
        
        print('Exist!')
        
    print('Gereing of fasta files ...')
    headers = []
    with open('{}/nrps_domains.fasta'.format(path), 'w') as fasta:
        for ind in table[table['Domain name'].str.contains('AMP-binding')].index:
            
            ID = table['ID'][ind]
            DOMAIN = table['Domain name'][ind]
            TRANSLATE = table['Sequence'][ind]        

            if '>{}_{}'.format(ID, DOMAIN) in headers:
                continue

            fasta.write('>{}_{}\n{}\n'.format(ID, 
                                            DOMAIN, 
                                            TRANSLATE))

            headers.append('>{}_{}'.format(ID, DOMAIN))

    print('Fasta files were generated successfully!')

    #Importtant to set on all hmms
    print('Searching values of sequense with HMM ...')
    hmms_s = listdir(hmms)

    for sub in hmms_s:
    
        call('hmmsearch -Z 1000 --cpu {} {}/{}/{}.hmm {}/nrps_domains.fasta > {}/{}.out'.format(cpu,
                                                                                                hmms,
                                                                                                sub,
                                                                                                sub,
                                                                                                path, 
                                                                                                hmm_out, 
                                                                                                sub), shell=True)

    print('Searching is done successfully!')
