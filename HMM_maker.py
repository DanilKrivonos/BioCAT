import argparse
from pandas import read_csv
from subprocess import call
from os import listdir, mkdir

def HMM_make(path, output, taxon, cpu):
    if taxon == 'bacteria':

        hmms = './HMM/Bacteria_HMM'

    else:

        hmms = './HMM/Fungi_HMM'

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

            if '>{}'.format(DOMAIN) in headers:
                continue

            fasta.write('>{}\n{}\n'.format(DOMAIN, 
                                           TRANSLATE))

            headers.append('>{}'.format(DOMAIN))

    print('Fasta files were generated successfully!')

    #Importtant to set on all hmms
    print('Searching values of sequense with HMM ...')
    hmms_s = listdir(hmms)

    for sub in hmms_s:

        substrate = sub.split('_')[1][: -4]
        call('hmmsearch -Z 1000 --cpu {} {}/{} {}/nrps_domains.fasta > {}/{}.out'.format(cpu,
                                                                                        hmms,
                                                                                        sub,
                                                                                        path, 
                                                                                        hmm_out, 
                                                                                        substrate), shell=True)

    print('Searching is done successfully!')
