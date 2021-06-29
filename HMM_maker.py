import argparse
from pandas import read_csv
from subprocess import call
from os import listdir

def HMM_make(path, output, hmms='./HMM/'):

    table = read_csv(path + '/table.tsv', sep='\t')

    #making work directory to fasta files
    substrates = []
    call('mkdir {}/substrates_fasta'.format(path), shell=True)
    print('Gereing of fasta files ...')

    for ind in table[table['Domain name'].str.contains('AMP-binding')].index:
        
        substrates.append(table['Single aa prediction'][ind])
        ID = table['ID'][ind]
        DOMAIN = table['Domain name'][ind]
        SUB = table['Single aa prediction'][ind]
        TRANSLATE = table['Sequence'][ind]

        with open('{}/substrates_fasta/{}_{}.fasta'.format(path, DOMAIN, SUB), 'w') as fasta:
            
            fasta.write('>{}_{}_{}\n{}'.format(ID, 
                                            DOMAIN, 
                                            SUB, 
                                            TRANSLATE))
    print('Fasta files were generated successfully!')

    #Importtant to set on all hmms
    print('Searching values of sequense with HMM ...')
    hmm_out = '{}/HMM_results/'.format(output)
    call('mkdir {}'.format(hmm_out), shell=True)
    A_domains = listdir('{}/substrates_fasta/'.format(path))
    hmms_s = listdir(hmms)

    for dom in A_domains:
        for sub in hmms_s:
            
            call('hmmsearch {}/{}/{}.hmm {}/substrates_fasta/{} > {}/{}vs{}.out'.format(hmms, 
                                                                                    sub, 
                                                                                    sub,
                                                                                    path, 
                                                                                    dom, 
                                                                                    hmm_out, 
                                                                                    sub, 
                                                                                    dom[: -6]), shell=True)

    print('Searching is done successfully!')
