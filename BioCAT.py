import argparse 
import os 
import sys
import numpy as np
from json import load
from subprocess import call 
from pandas import read_csv
from PSSM_maker import PSSM_make
from HMM_maker import HMM_make
from antiSMASH_parser import generate_table_from_antismash
from Make_NRP_structure import parse_rBAN
from Combinatorics import skipper, pssm, shuffle_matrix, make_combine

parser = argparse.ArgumentParser(description='Pipeline, which help to find biosynthesis gene clusters of NRP')
parser.add_argument('-name', 
                    type=str, 
                    help='Name of molecule or list of names', 
                    default="Unknown")
parser.add_argument('-smiles', 
                    type=str,
                    help='Chemical formula in smiles format', 
                    default=None)
parser.add_argument('-file_smiles',
                    type=str,
                    help='If you want to find a lot of substances you can give a file with smiles',
                    default=None)
parser.add_argument('-rBAN',
                    type=str,
                    help='Put here rBAN json',
                    default=None)
parser.add_argument('-antismash',
                    type=str,
                    help='Put here antismashs json',
                    default=None)                    
parser.add_argument('-genome', 
                    type=str, 
                    help='Fasta file with nucleotide sequence', 
                    default=None)
parser.add_argument('-NRPS_type_C', 
                    type=str, 
                    help='Expected NRPS type C', 
                    default='A')
parser.add_argument('-dif_strand', 
                    type=str, 
                    help='If your putative cluster can containd different strands genes', 
                    default=None)
parser.add_argument('-skip', 
                    type=int, 
                    help='Count of possible skippikng', 
                    default=0)
parser.add_argument('-out', 
                    type=str, 
                    help='Output directory, example: ./MY_path', 
                    default='./BioCAT_output')
parser.add_argument('-delta', 
                    type=str, 
                    help='Delta between PSSM and sequence of peptide', 
                    default=3)
args = parser.parse_args()

if args.genome == None and args.antismash == None:
    
    print('Error: Give a fasta or an antismash json!')
    sys.exit()
    
if args.genome != None:
    genome = args.genome
    fasta = open(genome)
    
output = args.out    
skip = args.skip
delta = int(args.delta)
NRPS_type = args.NRPS_type_C
dif_strand = args.dif_strand
#list with all substrates from antiSMASH
subtrate_stack = ['ser', 'thr', 'dhpg', 'hpg', 'gly', 'ala', 'val', 'leu', 'ile', 'abu', 'iva', 'pro', 'pip', 'asp', 'asn',
                  'glu', 'gln', 'aad', 'phe', 'trp', 'phg', 'tyr', 'bht', 'orn', 'lys', 'arg', 'cys', 'dhb', 'sal', 'nan']

#making working directiry
try:
    if args.out != './':
    
        os.mkdir(output)

except FileExistsError:

    print('The output directory already exists')

if args.smiles is not None and args.file_smiles is not None:
    
    print('Give a smiles or file if smi format!')
    import sys
    sys.exit()
    
if args.smiles is not None:

    smile_list = [args.smiles]
    ids = [args.name]

elif args.file_smiles is not None:
            
    smile_list = []
    ids = []

    with open(args.file_smiles) as smi:
        for line in smi:

            smile_list.append(line.split('\t')[1].replace('\n', ''))
            ids.append(line.split('\t')[0])

#call antiSMASH to take BGC
print('Finding biosynthesis gene clusters with antiSMASH ...\n')
if args.antismash is None:
    
    anti_out = output + 'antismash_result/'

    try:
        
        os.mkdir(anti_out)

    except FileExistsError:

        print('The output directory already exists')

    call('antismash {} --cb-general--output-dir {} --genefinding-tool prodigal'.format(genome, anti_out), shell=True)
    json_path = anti_out + ('.').join(os.path.split(genome)[1].split('.')[0: -1]) + '.json'
else:

    json_path = args.antismash
#parsing antiSMASH output json
generate_table_from_antismash(json_path, output, dif_strand)
#making of fasta files and hmmserching 
HMM_make(output, output)

def get_ids(outp):
    with open(outp) as json:

        js = load(json)
        
    try:

        name = js['id']

    except:

        name = js[0]['id']

    return  name

if args.rBAN is not None:
    
    smile_list = ['smi']
    ids = [get_ids(args.rBAN)]
 
with open('{}/Results.bed'.format(output), 'w') as bad_out:

    bad_out.write('Chromosome ID\tCoordinates of cluster\tStrand\tSubstance\tMiBiG ID\tPutative NRP sequence\th-score\n')

    for smi in range(len(smile_list)):
        if args.rBAN is None:
            
            smiles = '"' + smile_list[smi] + '"'
            idsmi = ids[smi].replace(' ', '_').replace('(', '').replace(')', '')
        
        #run rBUN
        print('Hydrolizing of substrate with rBAN ...')
        if args.rBAN is None:

            rBAN_path = output + '/{}_peptideGraph.json'.format(idsmi)
            call('java -jar rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {}/{}_ -discoveryMode'.format(idsmi, smiles, output, idsmi), shell=True)
            print('java -jar rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {}/{}_ -discoveryMode'.format(idsmi, smiles, output, idsmi[1: -1]))
            
        else:

            rBAN_path = args.rBAN

        new_EP = parse_rBAN(rBAN_path, NRPS_type)
        #if structure doesnt contain aminoacids
        if new_EP is None:

            print('Unparseable structure!')
            break 

        #Mking list of monomers
        print('Buildung amino graph')
        [print('Amino sequence of your substance: {}\n'.format(new_EP[seq])) for seq in new_EP]

        aminochain = make_combine(new_EP, subtrate_stack)
        #making PSSMs
        PSSM_make(search = output + '/HMM_results/', aminochain=aminochain, out = output, delta=delta)
        #Importing all PSSMs
        folder =  output + '/PSSM/'
        ITER = 100
        files = os.listdir(folder)
        #Trying to find some vsiants of biosynthesis
        if len(files) == 0:

            check = 0
            NRPS_type = 'C'

            while check != 1:

                print('Peptide sequence exceeds cluster landing attachment\nTry to check type C NRPS...')
                new_EP = parse_rBAN(rBAN_path, NRPS_type)
                [print('Amino sequence of your substance: {}\n'.format(new_EP[seq])) for seq in new_EP]
                aminochain = make_combine(new_EP, subtrate_stack)
                PSSM_make(search = output + '/HMM_results/', aminochain=aminochain, out = output, delta=delta)
                folder =  output + '/PSSM/'
                files = os.listdir(folder)

                if len(files) != 0:
                    check = 1
                    break

                if dif_strand == 'Have':

                    NRPS_type = 'A'

                #Dont split cluster
                if len(files) == 0 and dif_strand is None:

                    print('Trying to find putative cluster from different strands ...')
                    dif_strand = 'Have'
                    generate_table_from_antismash(json_path, output, dif_strand)
                
                if NRPS_type == 'A' and dif_strand == 'Have':

                    check = 1

        #Recording final output
        table = read_csv(output + '/table.tsv', sep='\t')
        
        if len(files) == 0:

            print('Organism have no putative cluster')

        for file in files:
            
            try:

                BGC_ID = file.split('.')[0].split('_A_')[1]

            except:

                continue

            for ind in table[table['ID'].str.contains(BGC_ID)].index:

                Name = table[table['ID'].str.contains(file.split('.')[0].split('_A_')[1])]['Name'][ind]
                print(file.split('.')[0].split('_A_')[1])
                Coord_cluster = table['Coordinates of cluster'][ind]
                strand = table['Gen strand'][ind]

                break

            if '_A_' not in file:
                continue


            BGC = read_csv(folder + file, sep='\t')
            
            if skip == 0:
                
                BGC = [BGC]
                
            else:
                
                BGC == skipper(BGC, skip)
                
            for matrix in BGC:
                
                EPs = make_combine(new_EP, subtrate_stack, matrix, delta)
                print(EPs, '-------------------')
                if EPs is None:
                    continue
                #Calculating TP score
                for v in EPs:
                    
                    shuffled_scores = []
                    
                    if len(matrix) == 1:
                        continue
                    for i in range(ITER):
                        
                        shuffled_scores.append(pssm(v, shuffle_matrix(matrix)))
                    
                    shuffled_scores = np.array(shuffled_scores)
                    TEST = pssm(v, matrix)
                    target_score = pssm(v, matrix)                
                    prob = len(shuffled_scores[shuffled_scores < target_score])/len(shuffled_scores)             
                    bad_out.write('{}\n'.format('\t'.join([Name,
                                                           Coord_cluster, 
                                                           strand,
                                                           ids[smi], 
                                                           BGC_ID,
                                                           '--'.join(v),
                                                           str(prob)])))
            
print('Job is done!')
