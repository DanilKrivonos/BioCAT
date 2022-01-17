import argparse 
import os 
import sys
import logging
import numpy as np
from json import load
from subprocess import call 
from pandas import read_csv, DataFrame
from BioCAT.src.PSSM_maker import PSSM_make
from BioCAT.src.HMM_maker import HMM_make
from BioCAT.src.antiSMASH_parser import generate_table_from_antismash
from BioCAT.src.Make_NRP_structure import parse_rBAN
from BioCAT.src.Combinatorics import skipper, get_score, shuffle_matrix, make_combine, multi_thread_shuffling, get_minim_aminochain
from BioCAT.src.Exploration_mode import exploration_mode
from BioCAT.src.Technical_functions import parse_smi_file, run_antiSMASH, run_rBAN, get_ids, check_pept_seq, get_results_to_csv, make_standard
from BioCAT.src.Calculating_scores import give_results

def main():
    parser = argparse.ArgumentParser(description='BioCAT is a tool, which estimates the' + 
                                    'likelihood that a given orgnism is capable of producing of a given NRP')
    group1 = parser.add_argument_group("Genome arguments")
    group2 = parser.add_argument_group("Chemical arguments")
    group3 = parser.add_argument_group("Technical arguments")
    group4 = parser.add_argument_group("Advanced arguments")
    # Group of genome arguments
    group1.add_argument('-antismash',
                        type=str,
                        help='antiSMASH *.json output file (either -antismash or -genome parameter should be specified)',
                        default=None)      
    group1.add_argument('-genome',
                        type=str,
                        help='Path to the fasta file with nucleotide sequence (either -antismash or -genome parameter should be specified)',
                        default=None)
    # Group of chemical arguments
    group2.add_argument('-name',
                        type=str,
                        help='Name of the given molecule (optional)',
                        default="Unknown")
    group2.add_argument('-smiles', 
                        type=str,
                        help='NRP chemical structure in the SMILES format (either -smiles or -file_smiles parameter should be specified)',
                        default=None)
    group2.add_argument('-file_smiles',
                        type=str,
                        help='.smi file with one or more NRPs. Each row should contain two ' + 
                        'columns: name of the NRP and SMILES string. Columns should be separated ' +  
                        'by tabulation. (Either -smiles or -file_smiles parameter should be specified.)',
                        default=None)
    group2.add_argument('-rBAN',
                        type=str,
                        help='rBAN peptideGraph.json output file',
                        default=None)
    group2.add_argument('-NRPS_type',
                        type=str,
                        help='Expected NRPS type (default A+B)',
                        default='A+B')
    # Group of technical arguments
    group3.add_argument('-iterations',
                        type=int,
                        help='Count of shuffling iterations (default 100)',
                        default=500)
    group3.add_argument('-delta',
                        type=str,
                        help='The maximum number of gaps in the molecule (default 3). Gaps are assigned as "nan".',
                        default=3)
    group3.add_argument('-threads',
                        type=int,
                        help='Number of threads (default 8)',
                        default=8)
    group3.add_argument('-out',
                        type=str,
                        help='Output directory (default ./BioCAT_output)',
                        default='./BioCAT_output')
    # Group of Advanced arguments
    group4.add_argument('-skip',
                        type=int,
                        help='Count of modules which can be skipped (default 0). ' + 
                        'Not recommended to use unless the user is sure about module skipping.',
                        default=0)
    group4.add_argument('--disable_pushing_type_B',
                        help='By defult, the algorithm tries to truncate peptide ' +
                            'fragments cutting edge monomers to find all possible identical ' +
                            'peptide fragments in the structure. If disabled, only ' +
                            'the identity of full peptide fragments is considered.',
                        action='store_true',
                        default=False)
    group4.add_argument('--disable_dif_strand',
                        help='By default, the protoclusters predicted by antiSMASH are ' +
                            'subdivided according to the assumption that each cluster ' +
                            'should contain only genes located on the same strand of the genome. ' +
                            'If disabled, protoclusters annotated by antiSMASH are ' +
                            'used as minimal clusters.',
                        action='store_true',
                        default=False)
    group4.add_argument('--disable_exploration',
                        help='By default, the algorithm tries to find the optimal ' +
                            'alignment combining alignment options in all possible ' +
                            'ways. If disabled, alignment is performed in according ' +
                            'to the given options strictly.',
                        action='store_true',
                        default=False)
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    # Parsing agruments
    cpu = args.threads
    disable_exploration = args.disable_exploration
    output = args.out    
    skip = args.skip
    delta = int(args.delta)
    NRPS_type = args.NRPS_type
    dont_dif_strand = args.disable_dif_strand
    genome = args.genome
    iterations = args.iterations
    off_push_type_B = args.disable_pushing_type_B
    if args.genome == None and args.antismash == None:
        
        print('Error: Give a fasta or an antismash json!')
        sys.exit()
    # Trying to make new directory
    try:
        if args.out != './':
        
            os.mkdir(output)

    except FileExistsError:

        print('The output directory already exists!')
    
    a_logger = logging.getLogger()
    a_logger.setLevel(logging.DEBUG)

    output_file_handler = logging.FileHandler("{}/BioCAT.log".format(output))
    stdout_handler = logging.StreamHandler(sys.stdout)

    a_logger.addHandler(output_file_handler)
    a_logger.addHandler(stdout_handler)
    # Checking have the programm have a smiles or file with smiles
    if args.smiles is not None and args.file_smiles is not None:
        
        a_logger.debug('Give a smiles or file if smi format!')
        sys.exit() # if have not aborting runs
    # Getting smiles and names of substances
    if args.rBAN == None:
        smile_list, ids = parse_smi_file(args.smiles, args.name, args.file_smiles)
    # Call antiSMASH to take BGC
    a_logger.debug('Finding biosynthesis gene clusters with antiSMASH ...\n')
    # Run antiSMASH if it neccecery or return path to directory 
    json_path = run_antiSMASH(args.antismash, output, genome, cpu)
    # Making antiSMASH dataes to table
    generate_table_from_antismash(json_path, output, dont_dif_strand)
    # Making of fasta files and hmmserching 
    HMM_make(output, output, cpu)
    # Getting substance name if we have only rBAN json
    a_logger.debug('Biosynthesis clusters have been successfully discovered!')
    if args.rBAN is not None:
        """
        If we have rBAN json we make empty smile_list,
        because we should not run rBAN yet and we can 
        only parse existing json output. Especialy, we
        gets name of the substance from peptideGraph.json.
        """
        smile_list = ['smi']
        ids = [get_ids(args.rBAN)]

    for smi in range(len(smile_list)):
        # Output info dictionary
        dict_res = {'Chromosome ID': [],
                'Coordinates of cluster': [],
                'Strand': [],
                'Substance': [],
                'BGC ID': [],
                'Putative linearized NRP sequence': [],
                'Biosynthesis profile': [],
                'Sln score': [],
                'Mln score': [],
                'Sdn score': [],
                'Mdn score': [],
                'Sdt score': [],
                'Mdt score': [],
                'Slt score': [],
                'Mlt score': [],
                'Relative score': [],
                'Binary': []
                }
        a_logger.debug('Calculating probability for {}'.format(ids[smi]))
        # Run rBAN if it neccecery or return path to directory 
        a_logger.debug('Buildung amino graph')
        rBAN_path = run_rBAN(args.rBAN, ids[smi], smile_list[smi], output)
        PeptideSeq = parse_rBAN(rBAN_path, NRPS_type, off_push_type_B)
        #if structure doesnt contain aminoacids
        if PeptideSeq is None:

            a_logger.debug('Unparseable structure!')
            get_results_to_csv(dict_res, output, ids[smi])# Make empty Results.csv
            continue 
        # Standartization of monomer names
        PeptideSeq = make_standard(PeptideSeq)  
        # Check correctness of the stcucture
        PeptideSeq = check_pept_seq(PeptideSeq)

        if PeptideSeq is None:

            a_logger.debug('Unparseable structure!')
            get_results_to_csv(dict_res, output, ids[smi])# Make empty Results.csv
            continue 
        
        # Calculating length of smaller variant
        length_min = get_minim_aminochain(PeptideSeq)
        for BS_type in PeptideSeq:

            a_logger.debug('Variants of fragments for Type {} biosynthetic path:'.format(BS_type))

            if PeptideSeq[BS_type] == {}:
    
                print('\tNone variants')
                continue

            for i in PeptideSeq[BS_type]:
                
                variant = str(PeptideSeq[BS_type][i]).replace('], [', ' & ')
                variant = variant.replace('[', '').replace(']', '').replace(', ', '--')
                variant = variant.replace("'", '')
                a_logger.debug('\t' + str(int(i) + 1) + ': ' + variant)
        
        # Making PSSMs
        PSSM_make(search = output + '/HMM_results/', aminochain=length_min, out = output, delta=delta, substance_name=ids[smi])
        #Importing all PSSMs
        folder =  output + '/PSSM_{}/'.format(ids[smi])
        files = os.listdir(folder)
        #Trying to find some vsiants of biosynthesis
        if disable_exploration == False:
            if len(files) == 0:

                PeptideSeq, NRPS_type = exploration_mode(rBAN_path, output, json_path, delta, substance_name=ids[smi])
                
        # Check availability of PSSMs
        files = os.listdir(folder)

        if len(files) == 0: 

            a_logger.debug('Organism has no putative cluster') 
            get_results_to_csv(dict_res, output, ids[smi])# Make empty Results.csv
            continue# If have no cluster -> brake it
        # Importing table with information about cluster
        table = read_csv(output + '/table.tsv', sep='\t')
        dict_res = give_results(dict_res, folder, files, table, ids[smi], PeptideSeq, skip, cpu, iterations)                                     
        #Recording Data Frame
        get_results_to_csv(dict_res, output, ids[smi])
        a_logger.debug('Results file for {} is recorded!'.format(ids[smi]))

    a_logger.debug('BioCAT processing is done!')
if __name__ == "__main__":
    main()
