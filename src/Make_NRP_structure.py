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

    tour_list_cop = tour_list.copy()

    for tour in tour_list:

        if tour[-1: : -1] in tour_list:
            tour_list.remove(tour[-1: : -1])

    return tour_list

def macthing_templates(templates_peptide, templates_minipept):

    templates_minipept_cop = templates_minipept.copy()

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

            templates_minipept_cop.remove(i)
        
    return templates_minipept_cop

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

    EP_cop = EP.copy()

    for var in EP_cop:

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
                
        if lost_edge[0] in numbers or lost_edge[1] in numbers:
            
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
# The function adding type B variant of biosynthesis to PeptideSeq
def compare_type_B_tour(add_tour, compare, PeptideSeq_cop, bios_path, key):
    if add_tour == compare:
        if add_tour not in PeptideSeq_cop[bios_path].values():
            if key not in PeptideSeq_cop['B']:

                PeptideSeq_cop['B'][key] = []

            if add_tour not in PeptideSeq_cop['B'][key]:
                
                PeptideSeq_cop['B'][key].append(add_tour) #because tour == tour x
            
    return PeptideSeq_cop
# Finding type B biosynthesis path
def type_B(PeptideSeq, Push='None'):
    
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

                    add_tour = PeptideSeq_cop[bios_path][var][tour]
                    compare = PeptideSeq_cop[bios_path][var][tourx]
                    new_PeptideSeq_copEP_cop = compare_type_B_tour(add_tour, compare, PeptideSeq_cop, bios_path, key)
                    # The push module need to approximate two homology NRP fragments, 
                    # which potentuly synthesys with type B variant
                    # Example BGC0000359: fuscachelin
                    if not Push is None: 
                        if len(PeptideSeq_cop[bios_path][var][tour]) != len(PeptideSeq_cop[bios_path][var][tourx]):
                            continue
                            
                        if PeptideSeq_cop[bios_path][var][tour] == PeptideSeq_cop[bios_path][var][tourx]:
                            continue
                        # In case if the last element will not match
                        add_tour = PeptideSeq_cop[bios_path][var][tour][: -1]

                        if add_tour == []:
                            continue
                    
                        compare = PeptideSeq_cop[bios_path][var][tourx][: -1]
                        PeptideSeq_cop = compare_type_B_tour(add_tour, compare, PeptideSeq_cop, bios_path, key)
                        # In case if the first element will not match
                        add_tour = PeptideSeq_cop[bios_path][var][tour][1: ]

                        if add_tour == []:
                            continue
                    
                        compare = PeptideSeq_cop[bios_path][var][tourx][1: ]
                        PeptideSeq_cop = compare_type_B_tour(add_tour, compare, PeptideSeq_cop, bios_path, key)
        key += 1

    return PeptideSeq_cop
    
def get_CPC(templates_minipept, js):
    # Make CPC(core peptide chain)
    CPC = [] 
    # Give it afunction 
    for tp in templates_minipept:
        for atom in js['atomicGraph']['atomicGraph']['atoms']:
            if atom['cdk_idx'] in tp:
                # N is conditional C end 
                # C is conditional N end 
                if atom['name'] == 'N':
                    
                    monomer_end = atom['matchIdx']
                    
                if atom['name'] == 'C':
                    
                    monomer_start = atom['matchIdx']
                    
        CPC.append([monomer_start, monomer_end]) #Directed edge fom conditional N to conditional C end
    return CPC

def get_peptide_elements(templates_peptide, js):
    
    core_peptide_elements = []
    
    for tp in templates_peptide:
        for atom in js['atomicGraph']['atomicGraph']['atoms']:
            if atom['cdk_idx'] in tp:
                if atom['matchIdx'] in core_peptide_elements:
                    continue
                    
                core_peptide_elements.append(atom['matchIdx']) #Adding aminoacids
    return core_peptide_elements

def get_non_classic(templates_non_classic, js):
    
    non_classic = []

    for tp in templates_non_classic:
        for atom in js['atomicGraph']['atomicGraph']['atoms']:
            if atom['cdk_idx'] in tp:
                if atom['matchIdx'] in non_classic:
                    continue

                non_classic.append(atom['matchIdx']) #Adding amino
    return non_classic

def get_AA(amino_acids_atoms, js):
    
    AAs = []

    for tp in amino_acids_atoms:
        for atom in js['atomicGraph']['atomicGraph']['atoms']:
            if atom['cdk_idx'] in tp:
                if atom['matchIdx'] in AAs:
                    continue

                AAs.append(atom['matchIdx']) #Adding amino
    return AAs

def get_direction(EP, CPC):
    
    for var in EP:
        for tour in EP[var]:
            if tour[-2 :] in CPC:
                continue

            EP[var][EP[var].index(tour)] = tour[: : -1]
            
    return EP

def add_acids(EP, amino_acids):

    EP_cop = EP.copy()

    for var in EP_cop:
        for tour in EP[var].copy():
            if len(tour) == 0:

                EP[var].remove(tour)
                continue
    
            for aa in amino_acids.copy():                    
                if aa not in tour:
                    continue

                amino_acids.remove(aa)

    if len(amino_acids) != 0:
        for var in EP:
            
            [EP[var].append([aa]) for aa in amino_acids]
            
    return EP
#Type C compression
def type_C(EP, js):
    
    repeats = Find_repeats(EP, js)
    EP_cop = EP.copy()

    for var in EP_cop:
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
# The function return list of monomer, which not AA, but pupular compounds
def find_not_aa_monomers(js):
    
    compare_dict = {'iva': 'CC(C)CC(=O)O',
                    'hiv': 'CC(C)C(C(=O)O)O',
                   # 'dhb': 'C1=CC(=C(C(=C1)O)O)C(=O)O', #dOH-Bz
                    'pip': 'C1CCNC(C1)C(=O)O' # Hpr
                   }
    not_aa_monomer = []
    
    for i in js['monomericGraph']['monomericGraph']['monomers']:

        target = Chem.MolFromSmiles(i['monomer']['monomer']['smiles'])
        
        if target is None: # In some cases can be unparseble structure in rBAN output
            continue

        for mon in compare_dict.keys():
                
            compare = Chem.MolFromSmiles(compare_dict[mon])
            Tanimoto = BulkTanimotoSimilarity(Chem.RDKFingerprint(target), [Chem.RDKFingerprint(compare)])[0]
            
            if Tanimoto == 1:
                if mon == 'hiv':
                    
                    mon = 'iva' # iva == hiv, the most common case is hiv, but we have only iva hmm (homology structure)
                    
                not_aa_monomer.append(i['monomer']['index'])
                i['monomer']['monomer']['monomer'] = mon
                
    return not_aa_monomer, js
def parse_rBAN(outp, NRPS_type, push_B=None):
    
    peptide_bond = Chem.MolFromSmiles('C(=O)CN(C(=O)CN)')
    mini_pept = Chem.MolFromSmiles('NC=O')
    N_OH = Chem.MolFromSmiles('C(=O)CN(O)(C(=O)CN)')
    alpha_amino = Chem.MolFromSmiles('C(=O)CN')
    N_N = Chem.MolFromSmiles('C(=O)CNC(=O)CNN')
    N_P = Chem.MolFromSmiles('C(=O)CN(P)(C(=O)CN)')
    N_Cl = Chem.MolFromSmiles('C(=O)CN(Cl)(C(=O)CN)')
    #Pattern length 
    pept_pattern_length = len(peptide_bond.GetAtoms())
    N_OH_length = len(N_OH.GetAtoms())
    
    with open(outp) as json:

        js = load(json)
        
    substance = Chem.MolFromSmiles(js['isomericSmiles'])
    templates_peptide = substance.GetSubstructMatches(peptide_bond)
    templates_non_classic = list(substance.GetSubstructMatches(N_OH))
    templates_non_classic += list(substance.GetSubstructMatches(N_N))
    templates_non_classic += list(substance.GetSubstructMatches(N_P))
    templates_non_classic += list(substance.GetSubstructMatches(N_Cl))
    templates_minipept = list(map(list, substance.GetSubstructMatches(mini_pept)))
    templates_minipept = macthing_templates(templates_peptide, templates_minipept)
    amino_acids_atoms = substance.GetSubstructMatches(alpha_amino)
    # Make CPC(core peptide chain)
    CPC = get_CPC(templates_minipept, js)
    # Getting main participant of core peptide chain
    core_peptide_elements = get_peptide_elements(templates_peptide, js)
    # Finding amino acids in the case of modificated peptide bonds
    amino_acids = get_AA(amino_acids_atoms, js)
    # Check non classic peptide bonds 
    non_pept = []
    
    if len(templates_non_classic) > 0:
        
        non_classic = get_non_classic(templates_non_classic, js)
        CPC_cop = CPC.copy()
        
        for edge in CPC_cop:
            if edge[0] in non_classic and edge[1] in non_classic:
                
                CPC.remove(edge)
                non_pept.append(edge)
        
    EP = {}
    EP[0] = []
    EP[0].extend(find_eulerian_tour(CPC.copy()))
    #for cyclic peptides bonded only with peptide bonds
    if EP[0] == []:

        EP = cyclic_peptide(CPC)

        if EP == {}:
            EP = {0: [[]]}
    #Check N-ends
    EP = get_direction(EP, CPC)

    if len(non_pept) > 0:

        EP = add_non_pept(EP, non_pept) #Return edges with non classic peptide elements
    # Adding not aa popular monomer 
    not_aa_monomer, js = find_not_aa_monomers(js)
    # Append standing separately aminoacids
    EP = add_acids(EP, amino_acids.copy())
    # If it unparsible analisys
    if len(EP[0][0]) == 0 or len(EP[0]) == 0:
        
        return None
    
    else:
        
        EP = find_amino_acid(EP, amino_acids.copy()) # Removing edges with non amino acids
        # Remove empty tours and append standing separately aminoacids
        EP = add_acids(EP, amino_acids.copy())
        EP = add_acids(EP, not_aa_monomer) # add non AA monomers
        # Giving biosynthetic way name
        PeptideSeq = {'A': EP, 'B': {}, 'C': {}}
        PeptideSeq['A'] = get_monomer_names(PeptideSeq['A'], js['monomericGraph']['monomericGraph']['monomers'])
        if 'C' in NRPS_type:

            PeptideSeq['C'] = type_C(EP, js)
            PeptideSeq['C'] = get_monomer_names(PeptideSeq['C'], js['monomericGraph']['monomericGraph']['monomers'])

        if 'B' in NRPS_type:
            
            PeptideSeq = type_B(PeptideSeq, push_B)

        return PeptideSeq