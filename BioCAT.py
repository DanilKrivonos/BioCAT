import argparse 
import os 
import random
import numpy as np
import sys
from rdkit import Chem
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
parser.add_argument('-rBAN',
                    type=str,
                    help='Put here rBAN json',
                    default="No")
parser.add_argument('-antismash',
                    type=str,
                    help='Put here antismashs json',
                    default="No")                    
parser.add_argument('-genome', 
                    type=str, 
                    help='Fasta file with nucleotide sequence', 
                    default=None)
parser.add_argument('-NRPS_type', 
                    type=str, 
                    help='Choose potentioanl NRPS type: A, B, C', 
                    default='A')
parser.add_argument('-skip', 
                    type=int, 
                    help='Count of possible skippikng', 
                    default=0)
parser.add_argument('-out', 
                    type=str, 
                    help='Output directory, example: ./MY_path', 
                    default='./')
parser.add_argument('-delta', 
                    type=str, 
                    help='Delta between PSSM and sequence of peptide', 
                    default=3)
args = parser.parse_args()

if args.genome == None and args.antismash == 'No':
    
    print('Error: Give a fasta or an antismash json!')
    sys.exit()
    
output = args.out + '/BioCAT_output/'

if args.genome != None:
    genome = args.genome
    fasta = open(genome)
    
skip = args.skip
NRPS_type = args.NRPS_type
#list with all substrates from antiSMASH
subtrate_stack = ['ser', 'thr', 'dhpg', 'hpg', 'gly', 'ala', 'val', 'leu', 'ile', 'abu', 'iva', 'pro', 'pip', 'asp', 'asn',
                  'glu', 'gln', 'aad', 'phe', 'trp', 'phg', 'tyr', 'bht', 'orn', 'lys', 'arg', 'cys', 'dhb', 'sal', 'nan']
#announce funtions
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
        repeat = N_cnt
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

    return sequences

def make_combine(sequences, pssm='None'):
    for var in sequences:

        sequences[var] = make_standard(sequences[var])

    N_refer = 0
    for var in sequences:
        for cont in sequences[var]:
            for mon in cont:
                if mon not in subtrate_stack:

                    cont[cont.index(mon)] = 'nan'
    print(sequences)
    for i in sequences[0]:
        for num in i:

            N_refer += 1 

    if type(pssm) == DataFrame:
        

        answ = []
        for var in sequences:
            
            nx = []
            
            for i in sequences[var]:
                
                nx.append(i)
            
            answ.extend(create_variants(nx, len(pssm)))
        answ = list(set(answ))
        return answ
        
    else:
        
        return  N_refer
    
#******************************************* GETTING LINER STRUCTURE **********************************
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
def macthing_templates(templates_peptide, templates_minipept):
    for i in templates_minipept:    
        
        check = 0
        
        for tp in templates_peptide:
            if check == len(i):
                continue

            check = 0

            for si in i:
                if si in tp:

                    check+=1

        if check != len(i):

            templates_minipept.remove(i)
        
    return templates_minipept

def build_grpah(js):
    
    edge_graph = []
    bonds_atoms = {}
    c_end_atoms = {}
    
    for i in js['monomericGraph']['monomericGraph']['bonds']:

        edge_graph.append(i['bond']['monomers'])

        for name in i['bond']['monomers']:
            if name not in bonds_atoms:

                bonds_atoms[name] = []
                c_end_atoms[name] = []

            atms = i['bond']['monomers'].copy()
            atms.remove(name)
            neibor = atms[0]

            for si in js['atomicGraph']['atomicGraph']['bonds']:
                if i['bond']['atomicIndexes'][0] == si['cdk_idx']:

                    c_end_atoms[name].extend(si['atoms'])

            for si in js['monomericGraph']['monomericGraph']['monomers']:
                if si['monomer']['index'] == neibor:

                    bonds_atoms[name].extend(si['monomer']['atoms'])
                    
    return edge_graph, c_end_atoms, bonds_atoms

def get_C_ends(js, c_end_atoms, templates, corboxy_pattern_lenth):
     
    C_ends = []
    pept_components = []
    atom_all = {}
    
    for line in js['monomericGraph']['monomericGraph']['monomers']:

        atoms = []
        atom_all[line['monomer']['index']] = []
        atom_all[line['monomer']['index']].extend(line['monomer']['atoms'])
        atoms.extend(line['monomer']['atoms'])
        atoms.extend(c_end_atoms[line['monomer']['index']])
        
        for tp in templates:
                
            s_atoms = 0

            for idx in tp:
                if idx in atoms:

                    s_atoms += 1
                    
            if s_atoms == corboxy_pattern_lenth:
                break
    
        try:                        
            if s_atoms == corboxy_pattern_lenth:

                C_ends.append(line['monomer']['index'])
                
        except NameError:
            continue
            
    return C_ends, atom_all

def get_N_ends(js, bonds_atoms, templates_peptide, corboxy_pattern_lenth, pept_pattern_lenth):
    
    pept_components = []
    
    for line in js['monomericGraph']['monomericGraph']['monomers']:

        atoms_for_pept = []
        atoms_for_pept.extend(line['monomer']['atoms'])
        atoms_for_pept.extend(bonds_atoms[line['monomer']['index']])
        atoms_for_pept = set(atoms_for_pept)
        amino_acids = []


        for tp in templates_peptide:

            p_atoms = 0
                    
            for idx in tp:
                if idx in atoms_for_pept:
                    
                    p_atoms += 1

            if p_atoms == pept_pattern_lenth:
                
                pept_components.append(line['monomer']['index'])
                break

    return pept_components

def get_AA(js, bonds_atoms, amino_acids_atoms, alpha_amino_length):
    
    amino_acids = []
    
    for line in js['monomericGraph']['monomericGraph']['monomers']:

        atoms_for_pept = []
        atoms_for_pept.extend(line['monomer']['atoms'])
       # atoms_for_pept.extend(bonds_atoms[line['monomer']['index']])
        atoms_for_pept = set(atoms_for_pept)
        

        for tp in amino_acids_atoms:

            p_atoms = 0
                    
            for idx in tp:
                if idx in atoms_for_pept:
                    
                    p_atoms += 1

            if p_atoms == alpha_amino_length:
                
                amino_acids.append(line['monomer']['index'])
                break
            
    return amino_acids

def find_amino_acid(EP, amono_acids):
    for var in EP:
        for tour in EP[var]:
            if len(tour) == 1 and tour[0] not in amono_acids:
                
                cp_var = EP[var].copy()
                cp_var.remove(tour)
                EP[var] = cp_var
                continue
                
            if tour[0] not in amono_acids:
                
                cp_tour = tour
                cp_tour.remove(tour[0])
                EP[var][EP[var].index(tour)] = cp_tour
                
            if tour[-1] not in amono_acids:
                
                cp_tour = tour
                cp_tour.remove(tour[-1])
                EP[var][EP[var].index(tour)] = cp_tour
                
    return EP

def cutter(edge_graph):
    variants = []

    for ind in range(len(edge_graph)):
        
        edge_gr_cop = edge_graph.copy()
        edge_gr_cop.remove(edge_gr_cop[ind])
        variants.append(edge_gr_cop)
        
    return variants

def get_peptide(edge_graph, pept_components, non_amino_graph):
    
    non_pept = []
    
    for edge in edge_graph:
        if edge[0] not in pept_components or edge[1] not in pept_components:
            #adding non peptide components
            non_pept.append(edge)
            
    #removeing non peptide chain components        
    for non in non_pept:
        
        edge_graph.remove(non)
    for edge in edge_graph:
        if edge in non_amino_graph:
            
            edge_graph.remove(edge)
            
    return edge_graph, non_pept

def cyclic_peptide(edge_graph):
    
    EP = {}
    v = cutter(edge_graph)
    idx = 0

    for t in v:

        EP[idx] = []
        EP[idx].extend(find_eulerian_tour(t))
        idx += 1
        
    return EP

def add_non_pept(EP, non_pept):
    
    for var in EP:
        
        check = 0
 
        for tour in EP[var]:
            if len(tour) == 1:
                continue
             
            for lost_edge in non_pept:
                if lost_edge[1] not in tour:
                    
                    EP[var].append([lost_edge[1]])
                    check = 1
                    
                elif lost_edge[0] not in tour:

                    EP[var].append([lost_edge[0]])
                    check = 1

        if check == 0:

            EP[var].append(lost_edge)
    return EP

def get_monomer_names(EP, space):
    new_EP = {}

    for var in EP:
        ind = 0
        new_EP[var] = []
        for tour in EP[var]:

            new_EP[var].append([])

            for i in tour:

                for mono in space:

                    sub = mono['monomer']['monomer']['monomer']
                    index = mono['monomer']['index']
            
                    if index == i:

                        new_EP[var][ind].append(sub.lower())
            ind += 1
            
    return new_EP
#Split product on one part of NRPS synthesis
def Type_B(new_EP):
    
    new_EP_cop = new_EP.copy()
    
    for var in new_EP:
        for tour in new_EP[var]:
            for tourx in new_EP[var]:
                if tourx == tour:
                    if [tour] not in new_EP.values():
                        
                        new_EP_cop[(len(new_EP))] = [tour] #because tour == tour x
                        
    return new_EP_cop
#check N-end atom 
def N_check(EP, tmp_names, atom_all):
    for var in EP:
        for tour in EP[var]:
            for N_number in tmp_names:
                if N_number in atom_all[tour[0]]:

                    EP[var][EP[var].index(tour)] = tour[: : -1]
                    break
    return EP
#check C-end atom 
def C_check(EP, C_ends):
    if C_ends != []:
            for var in EP:
                for tour in EP[var]:
                    if tour[0] in C_ends:

                        EP[var][EP[var].index(tour)] = tour[ : : -1]
    return EP

def parse_rBAN(outp, NRPS_type, subtrate_stack):
    
    corboxy_pattern = Chem.MolFromSmiles('C(N)C(=O)O')
    peptide_bond = Chem.MolFromSmiles('C(=O)CN(C(=O)CN)')
    mini_pept = Chem.MolFromSmiles('NC=O')
    alpha_amino = Chem.MolFromSmiles('C(=O)CN')
    #Pattern length 
    corboxy_pattern_lenth = len(corboxy_pattern.GetAtoms())
    pept_pattern_lenth = len(peptide_bond.GetAtoms())
    alpha_amino_length = len(alpha_amino.GetAtoms())
    with open(outp) as json:

        js = load(json)
        substance = Chem.MolFromSmiles(js['isomericSmiles'])
        templates = substance.GetSubstructMatches(corboxy_pattern)
        templates_peptide = substance.GetSubstructMatches(peptide_bond)
        tmp_names = []
        templates_minipept = list(map(list, substance.GetSubstructMatches(mini_pept)))
        amino_acids_atoms = list(map(list, substance.GetSubstructMatches(alpha_amino)))
        templates_minipept = macthing_templates(templates_peptide, templates_minipept)

        for tp in templates_minipept:
            for idx in js['atomicGraph']['atomicGraph']['atoms']:

                [tmp_names.append(i) for i in tp if idx['cdk_idx'] == i and idx['name'] == 'N']

        edge_graph, c_end_atoms, bonds_atoms = build_grpah(js)
        C_ends, atom_all = get_C_ends(js, c_end_atoms, templates, corboxy_pattern_lenth)
        pept_components = get_N_ends(js, bonds_atoms, templates_peptide, corboxy_pattern_lenth, pept_pattern_lenth)
        amino_acids = get_AA(js, bonds_atoms, amino_acids_atoms, alpha_amino_length)
        edge_graph = list(map(list, set(map(tuple, edge_graph))))
        edge_graph, non_pept = get_peptide(edge_graph, pept_components)
        EP = {}
        EP[0] = []
        EP[0].extend(find_eulerian_tour(edge_graph.copy()))

        #for cyclic peptides bonded only with peptide bonds
        if EP[0] == []:

            EP = cyclic_peptide(edge_graph)

        else:

            EP = C_check(N_check(EP, tmp_names, atom_all), C_ends)

        if len(non_pept) > 0:

            EP = add_non_pept(EP, non_pept)
        #find amino acids
        EP = find_amino_acid(EP, amino_acids)
        new_EP = get_monomer_names(EP, js['monomericGraph']['monomericGraph']['monomers'])
        new_EP = Type_B(new_EP)

    return new_EP
# ********************************************************************************************************
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

#Skip function 
def skipper(pssm, skip):
    skips_fragments = []
    skip_steps = list(permutations(pssm.index, skip))
    skip_steps = list(map(list, skip_steps))
    [com.sort() for com in skip_steps]
    skip_steps = list(map(list, list(set(map(tuple, skip_steps)))))

    for step in skip_steps:

        pssm_step = pssm
        if skip != 0:
            for skip_module in step:

                pssm_step = pssm_step[pssm_step.index != skip_module]
        skips_fragments.append(pssm_step)
        
    return skips_fragments
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

    bad_out.write('Name of organism\tPotentialy MiBiG ID\tSubstance\tCoordinates of cluster\tStrand\tPresumptive sequence\th-score\n')

    for smi in range(len(smile_list)):

        smiles = '"' + smile_list[smi] + '"'
        idsmi = '"' + ids[smi] + '"'
        #run rBUN
        print('Hydrolizing of substrate with rBAN ...')
        if args.rBAN == 'No':
            
            call('java -jar rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {} -discoveryMode'.format(idsmi, smiles, output), shell=True)
            new_EP = parse_rBAN(output + '/peptideGraph.json', NRPS_type, subtrate_stack)
            
        else:
            
            new_EP = parse_rBAN(args.rBAN, NRPS_type, subtrate_stack)
            
        #Mking list of monomers
        print('Buildung amino graph')
        print(output)
        
        [print('Amino sequence of your substance: {}\n'.format(new_EP[seq])) for seq in new_EP]
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
        
        if len(files) == 0:

            print('Organism have no potentianal cluster')

        for file in files:
            
            try:

                BGC_ID = file.split('.')[0].split('_A_')[1]

            except:

                continue

            for ind in table[table['ID'].str.contains(BGC_ID)].index:

                Name = table[table['ID'].str.contains(file.split('.')[0].split('_A_')[1])]['Name'][ind]
                Coord_cluster = table['Coordinates of cluster'][ind]
                strand = table['Strand'][ind]

                break

            if '_A_' not in file:
                continue


            BGC = read_csv(folder + file, sep='\t')
            
            if skip == 0:
                
                BGC = [BGC]
                
            else:
                
                BGC == skipper(BGC, skip)
                
            for matrix in BGC:
                
                EPs = make_combine(new_EP, matrix)
                print(EPs)
                #Calculating TP score

                

                for v in EPs:

                    shuffled_scores = []

                    for i in range(ITER):

                        shuffled_scores.append(pssm(v, shuffle_matrix(matrix)))

                    shuffled_scores = np.array(shuffled_scores)
                    print(v)
                    TEST = pssm(v, matrix)
                    target_score = pssm(v, matrix)                
                    prob = len(shuffled_scores[shuffled_scores < target_score])/len(shuffled_scores)
                    bad_out.write('{}\n'.format('\t'.join([BGC_ID,
                                                            Name,
                                                            ids[smi],
                                                            Coord_cluster, 
                                                            strand, 
                                                            '--'.join(v),
                                                            str(prob)])))
            
print('Job is done!')
