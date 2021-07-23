from os import listdir
from Make_NRP_structure import parse_rBAN
from Combinatorics import get_minim_aminochain, make_standard
from antiSMASH_parser import generate_table_from_antismash
from PSSM_maker import PSSM_make

# The exploratin mode try to find every possible combinations, with the exception of skipping
# 1 - Try to find type C sequence 
# 2 - Try don't split cluster on subcluster with type A+B+C biosynthetic path
# 2 - Try don't split cluster on subcluster with type A+B biosynthetic path

def exploration_mode(rBAN_path, output, json_path, delta, hmm=None):

    check = 0
    NRPS_type = 'A+B+C'
    dif_strand = None

    while check != 1:
        print('Peptide sequence exceeds cluster landing attachment\nTry to check type C NRPS...')
        PeptideSeq = parse_rBAN(rBAN_path, NRPS_type)

        for BS_type in PeptideSeq:

            print('For {} biosynthetic path'.format(BS_type))
            [print('Variant of amino sequence of your substance: {}\n'.format(PeptideSeq[BS_type][seq])) for seq in PeptideSeq[BS_type]]
###############################********************************TO TESTING**************************************************************
        # Standartization of monomer names
        PeptideSeq = make_standard(PeptideSeq)
        # Calculating length of smaller variant
        length_min = get_minim_aminochain(PeptideSeq)
        #PSSM_make(search = args.hmm + '/HMM_results/', aminochain=length_min, out = output, delta=delta)
        PSSM_make(search = output + '/HMM_results/', aminochain=length_min, out = output, delta=delta)

###############################********************************TO TESTING**************************************************************
        folder =  output + '/PSSM/'
        files = listdir(folder)

        if len(files) != 0:
            
            check = 1

        if dif_strand == 'Have':
            
            NRPS_type = 'A+B'
            dif_strand = None
            print('Trying to find putative cluster from different strands ...')
            generate_table_from_antismash(json_path, output, dif_strand)
        #Dont split cluster
        if len(files) == 0 and dif_strand is None:
            if NRPS_type != 'A+B':

                print('Trying to find putative cluster from different strands ...')
                dif_strand = 'Have'
                generate_table_from_antismash(json_path, output, dif_strand)
        print(NRPS_type)
        print(dif_strand)
        if NRPS_type == 'A+B' and dif_strand == None:
            check =1