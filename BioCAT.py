import argparse 
import os 
import random
import numpy as np
import sys
from random import shuffle
from subprocess import call 
from json import load
from pandas import read_csv, DataFrame
from itertools import permutations, product
from PSSM_maker import PSSM_make
from HMM_maker import HMM_make
from antiSMASH_parser import anti_parse

parser = argparse.ArgumentParser(description='Pipeline, which help to find biosynthesis gene clusters of NRP')
parser.add_argument('-name', 
                    type=str, 
                    help='Name of molecule or list of names', 
                    default="Unknown")
parser.add_argument('-smiles', 
                    type=str,
                    help='Chemical formula in smiles format', 
                    default="Unknown")
parser.add_argument('-file_smiles',
                    type=str,
                    help='If you want to find a lot of substances you can give a file with smiles(example \n>NAME OF SUBSTANCE\nSMILES_FORMULA)',
                    default="Unknown")
parser.add_argument('-antismash',
                    type=str,
                    help='Put here antismashs json',
                    default="No")                    
parser.add_argument('-genome', 
                    type=str, 
                    help='Fasta file with nucleotide sequence', 
                    default=None)
parser.add_argument('-out', 
                    type=str, 
                    help='Output directory, example: ./MY_path', 
                    default='./')
parser.add_argument('-delta', 
                    type=str, 
                    help='Output directory, example: ./MY_path', 
                    default=3)
args = parser.parse_args()

if args.genome == None and args.antismash == 'No':
    
    print('Error: Give a fasta or an antismash json!')
    sys.exit()

output = args.out + '/BioCAT_output/'

if args.genome != None:
    genome = args.genome
    fasta = open(genome)

#announce funtions
#Findind Euler tour
def iterator(current_node, graph, tour):
    priv_t = tour.copy()
    
    copy = graph.copy()
    for edge in copy:
        if current_node in edge:

            c_edge = edge.copy()
            c_edge.remove(current_node)
            current_node = c_edge[0]
            
            if current_node not in tour:
                
                tour.append(current_node)

    if priv_t != tour:
        
        iterator(current_node, graph, tour)
        
    return tour

def find_eulerian_tour(graph):
    
    start = []
    for edge in graph:
        for node in edge:
            if node not in start:
                start.append(node)
            else:
                start.remove(node)
                
    start.sort()
    tour_list = []

    for st in start:

        tour = []
        current_node = st
        tour.append(current_node)
        tour = iterator(current_node, graph, tour)
        tour_list.append(tour)
        
    for tour in tour_list:
        if tour[-1: : -1] in tour_list:
            tour_list.remove(tour[-1: : -1])

    return tour_list

#Making combinations of modules
def create_variants(original_seq, len_place,):
    
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
        repeat=N_cnt
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

    variants = set(variants)

    return variants

#Standartization of monomers
def make_standard(sequences):
    
    check_1 = 0
    check_2 = 0
    for seq in sequences:
        for sub_AS in subtrate_stack:
            for sub in seq:

                if sub_AS in sub:
                    
                    sequences[sequences.index(seq)][seq.index(sub)] =  sub_AS
                    
                if 'x' in sub:
                    
                    sequences[sequences.index(seq)][seq.index(sub)] =  'nan'

            if sub_AS in seq[0]:

                check_1 = 1                

            if sub_AS in seq[-1]:

                check_2 = 1
        
        if check_1 == 0:
        
            seq.remove(seq[0])

        if check_2 == 0:
                
            seq.remove(seq[-2])

    return sequences

def make_combine(sequences, pssm='None'):
    
    sequences = make_standard(sequences)
    EP_vars = {}
    ind = 0
    EP_vars[ind] = []
    ct = 0

    for seq in sequences:   

        new_seq = []
        
        for monom in seq:
            if monom in subtrate_stack:
                
                new_seq.append(monom)
                
            else:
                 
                new_seq = []
                
        ct += 1

        if new_seq[-1: : -1] not in EP_vars.values():
            if 'nan' in new_seq[0]:
                
                new_seq = new_seq[1: ]
                
            elif 'nan' in new_seq[-1]:
                
                new_seq = new_seq[: -1]
                
            EP_vars[ind].extend(new_seq)
            
        if ct != len(sequences):

            ind += 1
            EP_vars[ind] = []
    N_refer = 0
        
    for i in EP_vars:
        for num in EP_vars[i]:
            
            N_refer += 1 

    if type(pssm) == DataFrame:
        
        nx = []
        
        for i in EP_vars.keys():
        
            nx.append(EP_vars[i])
            
        return create_variants(nx, len(pssm))
    
        
    else:
        
        return  N_refer
def builder(n_end, concensus, est_graph, EP):
    
    while len(est_graph) != 0:
        for bond in est_graph:
            if n_end in bond:
                
                c_end = bond.copy()
                c_end.remove(n_end)
                c_end = c_end[0]
                est_graph.remove(bond)
                
        for contig in EP:
            if c_end == contig[0]:

                concensus.append(contig)
                c_end = contig[-1]
                EP.remove(contig)
                continue
            if c_end == contig[-1]:

                concensus.append(contig[-1: : -1])
                c_end = contig[-1: : -1][-1]
                EP.remove(contig)
                continue

        builder(c_end, concensus, est_graph, EP)

def parse_rBAN(outp):
    with open('{}/peptideGraph.json'.format(outp)) as json:

        js = load(json)
        edge_graph = []
        est_graph = []

        for i in js['monomericGraph']['monomericGraph']['bonds']:
            if 'AMINO' not in i['bond']['bondTypes']:
                
                est_graph.append(i['bond']['monomers'])
            
            if 'AMINO' not in i['bond']['bondTypes']:
                continue

            edge_graph.append(i['bond']['monomers'])
        
        EP = find_eulerian_tour(edge_graph.copy())
        est_list = []
        [est_list.extend(i) for i in est_graph]
        concensus = []
        #Finding true N-edn
        if len(EP) > 1:
            for i in EP:
                if i[-1] not in est_list:

                    n_end = i[0]
                    concensus.append(i[-1: : -1])
                    EP.remove(i)  
                    
                elif i[0] not in est_list:

                    n_end = i[-1]
                    concensus.append(i)
                    EP.remove(i)  
            
            builder(n_end, concensus, est_graph, EP)
        
        else:
            
            concensus = EP

        new_EP = []
        ind = 0
        for tour in concensus:
            
            new_EP.append([])
            
            for i in tour:
                
                for mono in js['monomericGraph']['monomericGraph']['monomers']:

                    sub = mono['monomer']['monomer']['monomer']
                    index = mono['monomer']['index']

                    if index == i:

                        new_EP[ind].append(sub.lower())
            ind += 1
    
    return new_EP
#Making shuffle funtions
def shuffle_matrix(pssm_profile):
    
    profile = pssm_profile.copy()
    cols = profile.columns[1:]
    
    for col in cols:
        
        col_vals = list(profile[col])
        shuffle(col_vals)
        profile[col] = col_vals
        
    return profile     
        
def pssm(seq, pssm_profile):
    
    
   
    TS = []
    cut = 0
    N = pssm_profile.shape[0]
    #Calculating C-score(cluster-score)
    
    while N - cut - len(seq) != -1:
        seq_cnt = 0
        target_sum = 0

        for line in pssm_profile.index:
            if line < cut:
                continue
            try:
                
                target_sum += pssm_profile[seq[seq_cnt]][line]
                
            except KeyError:
                
                
                target_sum += 0
                
            except IndexError:
                
                target_sum += 0

            seq_cnt += 1
        
        cut += 1
        TS.append(target_sum)
    
    return max(TS)
#list with all substrates from antiSMASH
subtrate_stack = ['ser', 'thr', 'dhpg', 'hpg', 'gly', 'ala', 'val', 'leu', 'ile', 'abu', 'iva', 'pro', 'pip', 'asp','asn', 
                  'glu', 'gln', 'aad', 'phe', 'trp', 'phg', 'tyr', 'bht', 'orn', 'lys', 'arg', 'cys', 'dhb', 'sal', 'nan']

#making working directiry
try:

    os.mkdir(output)

except FileExistsError:

    print('The output directory already exists')

if args.smiles != None:

    smile_list=[args.smiles]
    ids = [args.name]

else:

    smiles = args.smiles
    smile_list = []

    for line in smiles:

        [smile_list.append(smi) for smi in line.split('\n') if '>' not in smi]
        [ids.append(name) for name in line.split('\n') if '>' in name]

with open('{}/Results.bed'.format(output), 'w') as bad_out:

    bad_out.write('Name of organism\tPotentialy MiBiG ID\tSubstance\tCoordinates of cluster\tStrain\tPresumptive sequence\th-score\n')

    for smi in range(len(smile_list)):

        smiles = '"' + smile_list[smi] + '"'
        idsmi = '"' + ids[smi] + '"'
        #run rBUN
        print('Hydrolizing of substrate with rBAN ...')
        call('java -jar rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {} -discoveryMode'.format(idsmi, smiles, output), shell=True)

    

        #Mking list of monomers
        print('Buildung amino graph')

        new_EP = parse_rBAN(output)

        [print('Amino sequence of your substance: {}\n'.format('--'.join(seq))) for seq in new_EP]
        #call antiSMASH to take BGC
        print('Finding biosynthesis gene clusters with antiSMASH ...\n')
        if args.antismash == 'No':
            
            anti_out = output + 'antismash_result/'
            try:
                
                os.mkdir(anti_out)

            except FileExistsError:

                print('The output directory already exists')

            call('antismash {} --cb-general --cb-knownclusters --output-dir {} --genefinding-tool prodigal'.format(genome, anti_out), shell=True)
            json_path = anti_out + ('.').join(os.path.split(genome)[1].split('.')[0: -1]) + '.json'
        else:

            json_path = args.antismash
        #parsing antiSMASH output json
        anti_parse(json_path, output)
        #making of fasta files and hmmserching 
        HMM_make(output, output, hmms='./HMM/')
        aminochain = make_combine(new_EP)
        #making PSSMs
        PSSM_make(search = output + 'HMM_results/', aminochain=aminochain, out = output, delta=args.delta)
        #Importing all PSSMs
        folder =  output + 'PSSM/'
        ITER = 100
        files = os.listdir(folder)
        
        #Recording final output
        table = read_csv(output + 'table.tsv', sep='\t')

        for file in files:
            
            try:

                BGC_ID = file.split('.')[0].split('_A_')[1]

            except:

                continue

            for ind in table[table['ID'].str.contains(BGC_ID)].index:

                Name = table[table['ID'].str.contains(file.split('.')[0].split('_A_')[1])]['Name'][ind]
                Coord_cluster = table['Coordinates of cluster'][ind]
                strain = table['Strain'][ind]

                break

            if '_A_' not in file:
                continue

            try:

                BGC = read_csv(folder + file, sep='\t')
            
            except:
                continue

            EPs = make_combine(new_EP, BGC)
            print(EPs, '\\,[pkwp[vmpovoui')
            print(file, EPs, sep='\n')
            for v in EPs:
                
                TEST = pssm(v, BGC)

                #Calculating TP score
                target_score = pssm(v, BGC)
                shuffled_scores = []
             #   print('Calculating C-scores for {}'.format(file))

                for i in range(ITER):

                    shuffled_scores.append(pssm(v, shuffle_matrix(BGC)))

             #   print('Calculating h-scores for {}'.format(file))
                shuffled_scores = np.array(shuffled_scores)
                prob = len(shuffled_scores[shuffled_scores < target_score])/len(shuffled_scores)
                
                bad_out.write('{}\n'.format('\t'.join([BGC_ID,
                                                        Name,
                                                        ids[smi],
                                                        Coord_cluster, 
                                                        strain, 
                                                        '--'.join(v),
                                                        str(prob)])))
            
print('Job is done!')
