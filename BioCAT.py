import argparse 
import os 
import sys
import numpy as np
from json import load
from subprocess import call 
from pandas import read_csv, DataFrame
from PSSM_maker import PSSM_make
from HMM_maker import HMM_make
from antiSMASH_parser import generate_table_from_antismash
from Make_NRP_structure import parse_rBAN
from Combinatorics import skipper, get_score, shuffle_matrix, make_combine, multi_thread_shuffling, get_minim_aminochain, make_standard
from Exploration_mode import exploration_mode
from Technical_functions import parse_smi_file, run_antiSMASH, run_rBAN, get_ids
from Calculating_scores import give_results
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
parser.add_argument('-NRPS_type', 
                    type=str, 
                    help='Expected NRPS type^', 
                    default='A+B')
parser.add_argument('-dif_strand', 
                    type=str, 
                    help='If your putative cluster can containd different strands genes', 
                    default=None)
parser.add_argument('-exploration', 
                    type=bool, 
                    help='If you want to try every variants of biosynthesis', 
                    default=True)
parser.add_argument('-skip', 
                    type=int, 
                    help='Count of possible skippikng', 
                    default=0)
parser.add_argument('-out', 
                    type=str, 
                    help='Output directory, example: ./MY_path', 
                    default='./BioCAT_output')
parser.add_argument('-cpu', 
                    type=str, 
                    help='Multiple treadings', 
                    default=10)
###############################********************************TO TESTING**************************************************************
parser.add_argument('-hmm', 
                    type=str, 
                    help='Output directory, example: ./MY_path',
                    default=None)
###############################********************************TO TESTING**************************************************************
parser.add_argument('-delta', 
                    type=str, 
                    help='Delta between PSSM and sequence of peptide', 
                    default=3)
args = parser.parse_args()
# Parsing agruments 
cpu = args.cpu
exploration = args.exploration
output = args.out    
skip = args.skip
delta = int(args.delta)
NRPS_type = args.NRPS_type
dif_strand = args.dif_strand
genome = args.genome

if args.genome == None and args.antismash == None:
    
    print('Error: Give a fasta or an antismash json!')
    sys.exit()
 # Trying to make new directory
try:
    if args.out != './':
    
        os.mkdir(output)

except FileExistsError:

    print('The output directory already exists!')   
# Checking have the programm have a smiles or file with fmiles
if args.smiles is not None and args.file_smiles is not None:
    
    print('Give a smiles or file if smi format!')
    sys.exit() # if have not aborting runs
# Getting smiles and names of substances
if args.rBAN == None:
    smile_list, ids = parse_smi_file(args.smiles, args.name, args.file_smiles)
# Call antiSMASH to take BGC
print('Finding biosynthesis gene clusters with antiSMASH ...\n')
# Run antiSMASH if it neccecery or return path to directory 
json_path = run_antiSMASH(args.antismash, output, genome, cpu)
# Making antiSMASH dataes to table
generate_table_from_antismash(json_path, output, dif_strand)
# Making of fasta files and hmmserching 
###############################********************************TO TESTING**************************************************************
HMM_make(output, output, cpu)
###############################********************************TO TESTING**************************************************************
# Getting substance name if we have only rBAN json
if args.rBAN is not None:
    
    smile_list = ['smi']
    ids = [get_ids(args.rBAN)]
#Output info dictionary
bed_out = {'Chromosome ID': [],
            'Coordinates of cluster': [],
            'Strand': [],
            'Substance': [],
            'BGC ID': [],
            'Putative linearized NRP sequence': [],
            'Biosynthesis profile': [],
            'Sln score': [],
            'Mln score': [],
            'Slt score': [],
            'Mlt score': [],
            'Sdn score': [],
            'Mdn score': [],
            'Sdt score': [],
            'Mdt score': [],
            'Relative score': []
            }

for smi in range(len(smile_list)):
    # Run rBAN if it neccecery or return path to directory 
    rBAN_path = run_rBAN(args.rBAN, ids[smi], smile_list[smi], output)
    PeptideSeq = parse_rBAN(rBAN_path, NRPS_type)
    #if structure doesnt contain aminoacids
    if PeptideSeq is None:

        print('Unparseable structure!')
        break 

    #Mking list of monomers
    print('Buildung amino graph')
    for BS_type in PeptideSeq:

        print('For {} biosynthetic path'.format(BS_type))
        [print('Variant of amino sequence of your substance: {}\n'.format(PeptideSeq[BS_type][seq])) for seq in PeptideSeq[BS_type]]
    # Standartization of monomer names
    PeptideSeq = make_standard(PeptideSeq)
    # Calculating length of smaller variant
    length_min = get_minim_aminochain(PeptideSeq)
    #making PSSMs
###############################********************************TO TESTING**************************************************************
    # print( args.hmm + '/HMM_results/')
    #PSSM_make(search = args.hmm + '/HMM_results/', aminochain=length_min, out = output, delta=delta)
    PSSM_make(search = output + '/HMM_results/', aminochain=length_min, out = output, delta=delta)

###############################********************************TO TESTING**************************************************************
    #Importing all PSSMs
    folder =  output + '/PSSM/'
    files = os.listdir(folder)
    #Trying to find some vsiants of biosynthesis
    if exploration is True:
        if len(files) == 0:
###############################********************************TO TESTING**************************************************************
            #exploration(rBAN_path, output, json_path, args.hmm)
            exploration_mode(rBAN_path, output, json_path, delta)
###############################********************************TO TESTING**************************************************************
    # Check availability of PSSMs
    if len(files) == 0:

        print('Organism have no putative cluster') 
        break# If have no cluster -> brake it
    # Importing table with information about cluster
    table = read_csv(output + '/table.tsv', sep='\t')
    bed_out = give_results(bed_out, folder, files, table, ids, PeptideSeq, length_min, skip, cpu)                                     
#Recording Data Frame
bed_df = DataFrame(data=bed_out)
bed_df.to_csv('{}/Results.bed'.format(output), sep='\t')           
print('Job is done!')