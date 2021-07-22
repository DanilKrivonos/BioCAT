from json import load
from rdkit import Chem
from rdkit.DataStructs import BulkTanimotoSimilarity

#Findind Euler tour with modified Fleury's algorithm
def find_subgraph_tour(current_node, graph, tour):
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
        
        find_subgraph_tour(current_node, graph, tour)
        
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
        tour = find_subgraph_tour(current_node, graph, tour)
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

def build_graph(js):
    
    edge_graph = []
    non_amino_graph = []
    bonds_atoms = {}
    c_end_atoms = {}
    
    for i in js['monomericGraph']['monomericGraph']['bonds']:

        edge_graph.append(i['bond']['monomers'])
        
        if 'AMINO' not in i['bond']['bondTypes']:
            
            non_amino_graph.append(i['bond']['monomers'])
            
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
                    
    return edge_graph, c_end_atoms, bonds_atoms, non_amino_graph

def get_C_ends(js, c_end_atoms, templates_carboxy, corboxy_pattern_length):
     
    C_ends = []
    atom_all = {}
    
    for line in js['monomericGraph']['monomericGraph']['monomers']:

        atoms = []
        atom_all[line['monomer']['index']] = []
        atom_all[line['monomer']['index']].extend(line['monomer']['atoms'])
        atoms.extend(line['monomer']['atoms'])
        atoms.extend(c_end_atoms[line['monomer']['index']])
        
        for tp in templates_carboxy:
                
            s_atoms = 0

            for idx in tp:
                if idx in atoms:

                    s_atoms += 1
                    
            if s_atoms == corboxy_pattern_length:
                break
    
        try:                        
            if s_atoms == corboxy_pattern_length:

                C_ends.append(line['monomer']['index'])
                
        except NameError:
            continue
            
    return C_ends, atom_all

def get_N_OH(js, bonds_atoms, templates_N_OH, N_OH_length):
    
    N_OHs = []
    atom_all = {}
    
    for line in js['monomericGraph']['monomericGraph']['monomers']:

        atoms = []
        atom_all[line['monomer']['index']] = []
        atom_all[line['monomer']['index']].extend(line['monomer']['atoms'])
        atoms.extend(line['monomer']['atoms'])
        atoms.extend(bonds_atoms[line['monomer']['index']])
        
        for tp in templates_N_OH:
                
            s_atoms = 0

            for idx in tp:
                if idx in atoms:

                    s_atoms += 1
                    
            if s_atoms == N_OH_length:
                break
    
        try:                        
            if s_atoms == N_OH_length:

                N_OHs.append(line['monomer']['index'])
                
        except NameError:
            continue
            
    return N_OHs

def get_N_ends(js, bonds_atoms, templates_peptide, corboxy_pattern_length, pept_pattern_length):
    
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

            if p_atoms == pept_pattern_length:
                
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

def get_peptide(edge_graph, pept_components, non_amino_graph, C_ends):
    
    non_pept = []
    edge_graph_new = edge_graph.copy()
    
    for edge in edge_graph:
        if edge[0] not in pept_components or edge[1] not in pept_components:
            #adding non peptide components
            non_pept.append(edge)
            
    #removeing non peptide chain components        
    
    for non in non_pept:
        
        edge_graph_new.remove(non)
        
    edge_graph = edge_graph_new.copy()
    
    for edge in edge_graph:
        if edge in non_amino_graph:
            
            edge_graph_new.remove(edge)
            
        elif edge[0] in C_ends and edge[1] in C_ends:
            
            edge_graph_new.remove(edge)
            
    return edge_graph_new, non_pept

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
        numbers = []
        
        for tour in EP[var]:
            if len(tour) == 1 and tour != []:
                continue
                
            numbers += tour
            
        for lost_edge in non_pept:
            if lost_edge[1] not in numbers and [lost_edge[1]] not in EP[var]:
                
                EP[var].append([lost_edge[1]])
                check = 1

            elif lost_edge[0] not in numbers and [lost_edge[0]] not in EP[var]:

                EP[var].append([lost_edge[0]])
                check = 1
        
        if check == 0:

            EP[var].append(lost_edge)
        if [] in EP[var]:
            
            EP[var].remove([])
            
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
def Type_B(PeptideSeq):
    
    PeptideSeq_cop = PeptideSeq.copy()
    key = 0
    for bios_path in PeptideSeq_cop:
        if bios_path == 'B':
            continue

        for var in PeptideSeq_cop[bios_path]:
            for tour in range(len(PeptideSeq_cop[bios_path][var])): 

                indxs = list(range(len(PeptideSeq_cop[bios_path][var])))
                indxs.remove(tour)

                for tourx in indxs:
                    if PeptideSeq_cop[bios_path][var][tourx] == PeptideSeq_cop[bios_path][var][tour]:
                        if PeptideSeq_cop[bios_path][var][tour] not in PeptideSeq_cop[bios_path].values():
                            if key not in PeptideSeq_cop['B']:

                                PeptideSeq_cop['B'][key] = []

                            if PeptideSeq[bios_path][var][tour] in PeptideSeq_cop['B'][key]:
                                continue
                                
                            PeptideSeq_cop['B'][key].append(PeptideSeq[bios_path][var][tour]) #because tour == tour x
        key += 1

    return PeptideSeq_cop
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
#Type C compression
def type_C(EP, js):
    
    repeats = Find_repeats(EP, js)

    for var in EP:
        for tour in EP[var]:        
            if len(repeats[var][EP[var].index(tour)]) == 0:
                continue
                
            else:
                for rep in repeats[var][EP[var].index(tour)]:
                    for rep_node in rep[1: ]:
                        EP[var][EP[var].index(tour)].remove(rep_node)

    return EP
#finding reapits in aminochain tours
def Find_repeats(EP, js):
    
    reapets = {}
    
    for var in EP:
        
        reapets[var] = {}
        
        for tour in EP[var]:
            
            smiles_dict = {}

            for node in tour:
                for i in js['monomericGraph']['monomericGraph']['monomers']:
                    if i['monomer']['index'] == node:
                
                        smiles_dict[node] = Chem.MolFromSmiles(i['monomer']['monomer']['smiles'])
                        
            reapets[var][EP[var].index(tour)] = Get_compare(smiles_dict)
            


    return reapets
#finding Tanimotosimmularitis
def Get_compare(smiles_dict):
    
    reps = []
    mols = list(smiles_dict.keys())
    reps_count = 0
    
    for idx in range(len(mols) - 1):
        
        try:

            Tanimoto = BulkTanimotoSimilarity(Chem.RDKFingerprint(smiles_dict[mols[idx]]), [Chem.RDKFingerprint(smiles_dict[mols[idx + 1]])])[0]

        except:
            
            Tanimoto = 0

        if Tanimoto == 1:

            check = 0

            if len(reps) != 0:
                for rep in reps:
                    if mols[idx] in rep:

                        check = 1
                        
                if check == 1:

                    reps[reps_count].append(mols[idx + 1])

                else:

                    reps_count += 1
                    reps.append([])
                    reps[reps_count].extend([mols[idx], mols[idx + 1]])
            else:

                reps.append([mols[idx], mols[idx + 1]])                            

    return reps

def parse_rBAN(outp, NRPS_type):
    
    corboxy_pattern = Chem.MolFromSmiles('C(N)C(=O)O')
    peptide_bond = Chem.MolFromSmiles('C(=O)CN(C(=O)CN)')
    mini_pept = Chem.MolFromSmiles('NC=O')
    alpha_amino = Chem.MolFromSmiles('C(=O)CN')
    N_OH = Chem.MolFromSmiles('C(=O)CN(O)(C(=O)CN)')
    N_N = Chem.MolFromSmiles('C(=O)CNC(=O)CNN')
    N_P = Chem.MolFromSmiles('C(=O)CN(P)(C(=O)CN)')
    N_Cl = Chem.MolFromSmiles('C(=O)CN(Cl)(C(=O)CN)')
    #Pattern length 
    corboxy_pattern_length = len(corboxy_pattern.GetAtoms())
    pept_pattern_length = len(peptide_bond.GetAtoms())
    alpha_amino_length = len(alpha_amino.GetAtoms())
    N_OH_length = len(N_OH.GetAtoms())
    print(outp)
    with open(outp, 'r') as json:

        js = load(json)
        
    substance = Chem.MolFromSmiles(js['isomericSmiles'])
    templates_carboxy = substance.GetSubstructMatches(corboxy_pattern)
    templates_peptide = substance.GetSubstructMatches(peptide_bond)
    templates_N_OH = list(substance.GetSubstructMatches(N_OH))
    templates_N_OH += list(substance.GetSubstructMatches(N_N))
    templates_N_OH += list(substance.GetSubstructMatches(N_P))
    templates_N_OH += list(substance.GetSubstructMatches(N_Cl))
    tmp_names = []
    templates_minipept = list(map(list, substance.GetSubstructMatches(mini_pept)))
    amino_acids_atoms = list(map(list, substance.GetSubstructMatches(alpha_amino)))
    templates_minipept = macthing_templates(templates_peptide, templates_minipept)

    for tp in templates_minipept:
        for idx in js['atomicGraph']['atomicGraph']['atoms']:

            [tmp_names.append(i) for i in tp if idx['cdk_idx'] == i and idx['name'] == 'N']

    edge_graph, c_end_atoms, bonds_atoms, non_amino_graph = build_graph(js)
    C_ends, atom_all = get_C_ends(js, c_end_atoms, templates_carboxy, corboxy_pattern_length)
    pept_components = get_N_ends(js, bonds_atoms, templates_peptide, corboxy_pattern_length, pept_pattern_length)
    amino_acids = get_AA(js, bonds_atoms, amino_acids_atoms, alpha_amino_length)
    edge_graph, non_pept = get_peptide(edge_graph, pept_components, non_amino_graph, C_ends)

    if len(templates_N_OH) > 0:
        
        N_OH = get_N_OH(js, bonds_atoms, templates_N_OH, N_OH_length)
        for edge in edge_graph:
            if edge[0] in N_OH and edge[1] in N_OH:
                
                edge_graph.remove(edge)
                non_pept.append(edge)
        
    EP = {}
    EP[0] = []
    EP[0].extend(find_eulerian_tour(edge_graph.copy()))
    #for cyclic peptides bonded only with peptide bonds

    if EP[0] == []:

        EP = cyclic_peptide(edge_graph)
        
        if EP == {}:
            EP = {0: [[]]}
    else:

        EP = C_check(N_check(EP, tmp_names, atom_all), C_ends)

    if len(non_pept) > 0:

        EP = add_non_pept(EP, non_pept)
    #find amino acids
    if len(EP[0][0]) == 0:
    
        return None

    else:
        
        EP = find_amino_acid(EP, amino_acids)
        # Giveng biosynthetic way name
        PeptideSeq = {'A': EP, 'B': {}, 'C': {}}

        if 'C' in NRPS_type:

            PeptideSeq['C'] = type_C(EP, js)
            PeptideSeq['C'] = get_monomer_names(PeptideSeq['C'], js['monomericGraph']['monomericGraph']['monomers'])

        PeptideSeq['A'] = get_monomer_names(PeptideSeq['A'], js['monomericGraph']['monomericGraph']['monomers'])

        if 'B' in NRPS_type:
            
            PeptideSeq = Type_B(PeptideSeq)

        return PeptideSeq

