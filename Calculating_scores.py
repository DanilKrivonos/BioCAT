from numpy import array
from pickle import load
from pandas import read_csv
from Combinatorics import multi_thread_shuffling, multi_thread_calculating_scores, make_combine, get_score
# Importing random forest model
Rf = load(open('./model.dump', 'rb'))
# The function generate list of shuflled matrix
def make_shuffle_matrix(matrix, cpu):

    module_shuffling_matrix = multi_thread_shuffling(matrix, ShufflingType='module', iterations=100, threads=cpu)
    substrate_shuffling_matrix = multi_thread_shuffling(matrix, ShufflingType='substrate', iterations=100, threads=cpu)
    return module_shuffling_matrix, substrate_shuffling_matrix

# The fujnction finds suquence with maximum possible value, results from alignment
def get_MaxSeq(matrix, variant_seq):

    MaxSeq = []
    subs = matrix.keys()[1: ]
    # Find sequence, wich have maximum alignment score 
    for idx in matrix.index:

        MAX_value = max(list(matrix.iloc[idx][1:]))

        for key in subs:
            if matrix[key][idx] == MAX_value:

                MaxSeq.append(key) # If two smonomer have same value
                break
    # Making two variants of MaxSeq
    MaxSeq_full = MaxSeq.copy()
    MaxSeq_nan = MaxSeq.copy()
    
    for max_sub_idx in range(len(MaxSeq)):
        if variant_seq[max_sub_idx] == 'nan':

            MaxSeq_nan[max_sub_idx] = 'nan' # Adding nan to MaxSeq

    return MaxSeq_full, MaxSeq_nan
# The function gives an information about clusters
def get_cluster_info(table, BGC_ID, target_file):
    for ind in table[table['ID'].str.contains(BGC_ID)].index:

        Name = table[table['ID'].str.contains(target_file.split('.')[0].split('_A_')[1])]['Name'][ind]
        Coord_cluster = table['Coordinates of cluster'][ind]
        strand = table['Gen strand'][ind]
        break

    return Name, Coord_cluster, strand

# Calculate scores 
def calculate_scores(variant_seq, matrix, substrate_shuffling_matrix, module_shuffling_matrix, cpu):    
    # Finding suquence with maximum possible value, results from alignment
    MaxSeq_full, MaxSeq_nan = get_MaxSeq(matrix, variant_seq)
    # Calculating shuffled scores
    Sln_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_nan, substrate_shuffling_matrix, type_value='log', iterations=100, threads=cpu))
    Mln_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_nan, module_shuffling_matrix, type_value='log', iterations=100, threads=cpu))
    Slt_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_full, substrate_shuffling_matrix, type_value='log', iterations=100, threads=cpu))
    Mlt_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_full, module_shuffling_matrix, type_value='log', iterations=100, threads=cpu))
    Sdn_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_nan, substrate_shuffling_matrix, type_value=None, iterations=100, threads=cpu))
    Mdn_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_nan, module_shuffling_matrix, type_value=None, iterations=100, threads=cpu))
    Sdt_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_full, substrate_shuffling_matrix, type_value=None, iterations=100, threads=cpu))
    Mdt_shuffled_score = array(multi_thread_calculating_scores(MaxSeq_full, module_shuffling_matrix, type_value=None, iterations=100, threads=cpu))
    # Calculating scores for target sequence
    log_target_score = get_score(variant_seq, matrix, type_value='log')
    non_log_target_score = get_score(variant_seq, matrix, type_value=None)
    # Calculating features scores
    Sln_score = len(Sln_shuffled_score[Sln_shuffled_score < log_target_score])/len(Sln_shuffled_score)
    Mln_score = len(Mln_shuffled_score[Mln_shuffled_score < log_target_score])/len(Mln_shuffled_score)
    Slt_score = len(Slt_shuffled_score[Slt_shuffled_score < log_target_score])/len(Slt_shuffled_score)
    Mlt_score = len(Mlt_shuffled_score[Mlt_shuffled_score < log_target_score])/len(Mlt_shuffled_score)
    Sdn_score = len(Sdn_shuffled_score[Sdn_shuffled_score < non_log_target_score])/len(Sdn_shuffled_score)
    Mdn_score = len(Mdn_shuffled_score[Mdn_shuffled_score < non_log_target_score])/len(Mdn_shuffled_score)
    Sdt_score = len(Sdt_shuffled_score[Sdt_shuffled_score < non_log_target_score])/len(Sdt_shuffled_score)
    Mdt_score = len(Mdt_shuffled_score[Mdt_shuffled_score < non_log_target_score])/len(Mdt_shuffled_score)
    # Calculating Relative score
    
    Relative_score = Rf.predict_proba([[Sln_score, Mln_score, 
                                Slt_score, Mlt_score, 
                                Sdn_score, Mdn_score,
                                Sdt_score, Mdt_score]])[0][1]
    
    return Sln_score, Mln_score, Slt_score, Mlt_score, Sdn_score, Mdn_score, Sdt_score, Mdt_score, Relative_score

def give_results(bed_out, folder, files, table, ids, PeptideSeq, length_min, skip, cpu):
    
    for target_file in files:
            
        try:

            BGC_ID = target_file.split('.')[0].split('_A_')[1]

        except:

            continue

        if '_A_' not in target_file:
            continue

        Name, Coord_cluster, strand = get_cluster_info(table, BGC_ID, target_file) # Getting information about cluster
        BGC = read_csv(folder + target_file, sep='\t')
        # Skipping mode
        if skip == 0:
                
            BGC = [BGC]
                
        else:
                
            BGC == skipper(BGC, skip)

        for matrix in BGC:
            if len(matrix) == 1:
                continue

            for BS_type in PeptideSeq:# For every biosynthesis profile pathways

                EPs = make_combine(PeptideSeq[BS_type], length_min, matrix, delta=3)

                if EPs is None: # If length sequnce can't be scaled to cluster size
                    continue
                # Generating shuffling matrix
                module_shuffling_matrix, substrate_shuffling_matrix =  make_shuffle_matrix(matrix, cpu)
                
                for variant_seq in EPs:

                    Sln_score, Mln_score, Slt_score, Mlt_score, Sdn_score, Mdn_score, Sdt_score, Mdt_score, Relative_score = calculate_scores(variant_seq, matrix, substrate_shuffling_matrix, module_shuffling_matrix, cpu)
                    #Recordind dictionary 
                    bed_out['Chromosome ID'].append(Name)
                    bed_out['Coordinates of cluster'].append(Coord_cluster)
                    bed_out['Strand'].append(strand)
                    bed_out['Substance'].append(ids)
                    bed_out['BGC ID'].append(BGC_ID)
                    bed_out['Putative linearized NRP sequence'].append('--'.join(variant_seq))
                    bed_out['Biosynthesis profile'].append('Type {}'.format(BS_type))
                    bed_out['Sln score'].append(Sln_score) #shaffling substrates in matrix with log score and nan in maximally possible sequence
                    bed_out['Mln score'].append(Mln_score) #shaffling modules matrix with log score and nan in maximally possible sequence
                    bed_out['Slt score'].append(Slt_score) #shaffling substrates matrix with log score in maximally possible sequence
                    bed_out['Mlt score'].append(Mlt_score) #shaffling modules matrix with log score in maximally possible sequence
                    bed_out['Sdn score'].append(Sdn_score) #shaffling substrates matrix without log score and nan in maximally possible sequence
                    bed_out['Mdn score'].append(Mdn_score) #shaffling modules matrix without log score and nan in maximally possible sequence
                    bed_out['Sdt score'].append(Sdt_score) #shaffling substrates matrix without log score in maximally possible sequence
                    bed_out['Mdt score'].append(Mdt_score) #shaffling modules matrix without log score in maximally possible sequence
                    bed_out['Relative score'].append(Relative_score) #Final score
    return bed_out