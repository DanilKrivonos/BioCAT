from Make_NRP_structure import parse_rBAN
from Combinatorics import get_minim_aminochain, make_standard
from PSSM_maker import PSSM_make

# The exploratin mode try to find every possible combinations, with the exception of skipping
# 1 - Try to find type C sequence 
# 2 - Try don't split cluster on subcluster with type A+B+C biosynthetic path
# 2 - Try don't split cluster on subcluster with type A+B biosynthetic path

def exploration(rBAN_path, output, json_path, hmm=None):
    if exploration is True:
        if len(files) == 0:

            check = 0
            NRPS_type = 'A+B+C'

            while check != 1:

                print('Peptide sequence exceeds cluster landing attachment\nTry to check type C NRPS...')
                new_EP = parse_rBAN(rBAN_path, NRPS_type)
                [print('Amino sequence of your substance: {}\n'.format(new_EP[seq])) for seq in new_EP]
                new_EP = make_standard(new_EP)
                lenth_min = get_minim_aminochain(new_EP)
    ###############################********************************TO TESTING**************************************************************
                output
                #PSSM_make(search = args.hmm + '/HMM_results/', aminochain=lenth_min, out = output, delta=delta)
                PSSM_make(search = output + '/HMM_results/', aminochain=lenth_min, out = output, delta=delta)

    ###############################********************************TO TESTING**************************************************************
                folder =  output + '/PSSM/'
                files = os.listdir(folder)

                if len(files) != 0:
                    check = 1
                    break

                if dif_strand == 'Have':

                    NRPS_type = 'A+B'

                #Dont split cluster
                if len(files) == 0 and dif_strand is None:

                    print('Trying to find putative cluster from different strands ...')
                    dif_strand = 'Have'
                    generate_table_from_antismash(json_path, output, dif_strand)
                
                if NRPS_type == 'A+B' and dif_strand == 'Have':

                    check = 1