from os import listdir
from BioCAT.src.Make_NRP_structure import parse_rBAN
from BioCAT.src.Combinatorics import get_max_aminochain, get_minim_aminochain
from BioCAT.src.Technical_functions import make_standard
from BioCAT.src.antiSMASH_parser import generate_table_from_antismash
from BioCAT.src.PSSM_maker import PSSM_make

# The exploratin mode try to find every possible combinations, with the exception of skipping

def exploration_mode(rBAN_path, output, json_path, delta, substance_name):
    """
    Exploration mode is the function of searching of every commonly meating 
    variants of NRPS biosynthesis. It's try to find optimal way to find potantial
    cluster.
    1 - Try to find type A+B+C sequences with subclusters and push type B
    2 - Try don't split cluster on subcluster with type A+B+C biosynthetic path and push type B
    2 - Try don't split cluster on subcluster with type A+B biosynthetic path without push type B

    Parameters
    ----------
    rBAN_path : str
        Path to rBANs peptideGraph.json output.
    output : str
        Path to BioCAT output directory.
    json_path : str
        Path to antismash json output file.
    delta : int
        Differense between length of core peptide chain and size of cluster.
    substance_name : str
         Name of analyzing substance.
    Returns
    -------
    PeptideSeq : dict
        Changed PeptideSeq.
    NRPS_type : str
        Chnged NRPS type.
    """
    check = 0
    NRPS_type = 'A+B+C'
    off_push_B = False
    dont_dif_strand = False

    while check != 1:
        print('Peptide sequence exceeds cluster landing attachment\nTry to check type C NRPS...')
        PeptideSeq = parse_rBAN(rBAN_path, NRPS_type, off_push_B)
        
        for BS_type in PeptideSeq:

            print('For {} biosynthetic paths'.format(BS_type))
            print('Variant of amino sequence of your substance:', BS_type, PeptideSeq[BS_type])

        # Standartization of monomer names
        PeptideSeq = make_standard(PeptideSeq)
        # Calculating length of smaller variant
        length_min = get_minim_aminochain(PeptideSeq)
        PSSM_make(search = output + '/HMM_results/', aminochain=length_min, out = output, delta=delta, substance_name=substance_name)
        folder =  output + '/PSSM_{}/'.format(substance_name)
        files = listdir(folder)

        if len(files) != 0:
            
            check = 1
            break
        if dont_dif_strand == True:

            NRPS_type = 'A+B'
            off_push_B = True
            dont_dif_strand = False
            print('Trying to find putative cluster from different strands ...')
            generate_table_from_antismash(json_path, output, dont_dif_strand)

        files = listdir(folder)
        # Dont split cluster on subclusters
        if len(files) == 0 and dont_dif_strand == False:
            if NRPS_type != 'A+B':

                print('Trying to find putative cluster from different strands ...')
                dont_dif_strand = True
                generate_table_from_antismash(json_path, output, dont_dif_strand)
        
        if NRPS_type == 'A+B' and dont_dif_strand == False:

            check = 1
            
    return PeptideSeq, NRPS_type