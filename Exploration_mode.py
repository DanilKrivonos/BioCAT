from os import listdir
from Make_NRP_structure import parse_rBAN
from Combinatorics import make_standard, get_max_aminochain, get_minim_aminochain
from antiSMASH_parser import generate_table_from_antismash
from PSSM_maker import PSSM_make

# The exploratin mode try to find every possible combinations, with the exception of skipping
# 1 - Try to find type A+B+C sequences with subclusters and push type B
# 2 - Try don't split cluster on subcluster with type A+B+C biosynthetic path and push type B
# 2 - Try don't split cluster on subcluster with type A+B biosynthetic path without push type B

def exploration_mode(rBAN_path, output, json_path, delta, substance_name, taxon):

    check = 0
    NRPS_type = 'A+B+C'
    push_B = 'Push'
    dif_strand = None

    while check != 1:
        
        print('Peptide sequence exceeds cluster landing attachment\nTry to check type C NRPS...')
        PeptideSeq = parse_rBAN(rBAN_path, NRPS_type, push_B)
        
        for BS_type in PeptideSeq:

            print('For {} biosynthetic paths'.format(BS_type))
            print('Variant of amino sequence of your substance:', BS_type, PeptideSeq[BS_type])

        # Standartization of monomer names
        PeptideSeq = make_standard(PeptideSeq, taxon)
        # Calculating length of smaller variant
        length_min = get_minim_aminochain(PeptideSeq)

        PSSM_make(search = output + '/HMM_results/', aminochain=length_min, out = output, delta=delta, substance_name=substance_name, taxon=taxon)
        folder =  output + '/PSSM_{}/'.format(substance_name)
        files = listdir(folder)

        if len(files) != 0:
            
            check = 1
            break
        if dif_strand == 'Have':
            
            NRPS_type = 'A+B'
            push_B = None
            dif_strand = None
            print('Trying to find putative cluster from different strands ...')
            generate_table_from_antismash(json_path, output, dif_strand)

        files = listdir(folder)
        # Dont split cluster on subclusters
        if len(files) == 0 and dif_strand is None:
            if NRPS_type != 'A+B':

                print('Trying to find putative cluster from different strands ...')
                dif_strand = 'Have'
                generate_table_from_antismash(json_path, output, dif_strand)
        
        if NRPS_type == 'A+B' and dif_strand == None:

            check = 1
            
    return PeptideSeq, NRPS_type