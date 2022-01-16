import os
import logging
from os import listdir
from json import load
from subprocess import call
from pandas import DataFrame

a_logger = logging.getLogger()
a_logger.setLevel(logging.DEBUG)
"""
This if the technical function, which using to run external softs and getting metha information about 
anlizyng substance.
"""
#Getting smiles formula and names of substances
def parse_smi_file(smiles, name, file_smiles):
    """
    The function takes smiles formula or file in smi format
    and converts it into two lists with smiles and names of a
    substances.

    Parameters
    ----------
    smiles : str
        Chemical formula of substance in SMILES format.
    Returns
    -------
    smile_list : list
        List with smiles for every analyzing structures.
    ids : list
        List with names of analyzing structures.
    """
    smile_list = []
    ids = []
    
    if smiles is not None:
        
        smile_list = [smiles]
        ids = [name.replace(' ', '_').replace('(', '').replace(')', '')]
    
    elif file_smiles is not None:
        with open(file_smiles) as smi:
            for line in smi:
                #Replacing hard parsible elements. For example actinomycin D will actinomycin_D.
                smile_list.append(line.split('\t')[1].replace('\n', ''))
                ids.append(line.split('\t')[0].replace(' ', '_').replace('(', '').replace(')', ''))
    
    return smile_list, ids

# Run antiSMASH if it neccecery or return path to directory 
def run_antiSMASH(antismash, output, genome, cpu):
    """
    Returns the path to antiSMASH json otput file. If there are no any antismash results
    in the output folder, the function runs the Antismash procesing.

    Parameters
    ----------
    antismash : str or None
        Path to antiSMASH json. None, if user has given a genome.
    output : str
        Path to BioCAT output directory.
    genome : str or None
        Path to genome of popential producer.
    cpu : int
        Number of threads used.
    Returns
    -------
    json_path : list
        Path to antismash json output file.
    """
    if antismash is None: # If we have not existing json
        
        anti_out = output + '/antismash_result/'
        
        try:
            #Making directiry with antismash results. If it necessary.    
            os.mkdir(anti_out)
        
        except FileExistsError:

            a_logger.debug('The output directory already exists')
        # Run antiSMASH with loose mode
        a_logger.debug('Running antiSMASH ...')
        call('antismash {} --cb-general --output-dir {} --genefinding-tool prodigal --skip-zip-file --enable-nrps-pks --minimal --cpus {}'.format(genome, anti_out, cpu), shell=True)
        json_path = anti_out + ('.').join(os.path.split(genome)[1].split('.')[0: -1]) + '.json'

    else:

        json_path = antismash

    return json_path

# Getting substance name if we have only rBAN json
def get_ids(outp):
    """
    Return name of processed substance.

    Parameters
    ----------
    outp : str
        Path to rBANs peptideGraph.json output.
    Returns
    -------
    name : str
        Name of analyzing substance.
    """
    with open(outp, 'r') as json:

        js = load(json)
        
    try:

        name = js['id']

    except:

        name = js[0]['id']

    return  name
    
# Run rBAN if it neccecery or return path to directory 
def run_rBAN(rBAN, ID, SMI, output):
    """
    Runs rBAN, if it necessary, and returns the path to
    peptideGraph.json output.

    Parameters
    ----------
    rBAN : str
        Path to rBANs peptideGraph.json output.
    ID : str
        Name of analyzing substance.
    SMI : str
        SMILES of analyzing substance.
    output : str
        Path to BioCAT output directory.
    Returns
    -------
    rBAN_path : str
        Path to rBANs peptideGraph.json output file.
    """
    a_logger.debug('Hydrolizing of substrate with rBAN ...')
    if rBAN is None:

        smiles = '"' + SMI + '"'
        #idsmi = ID.replace(' ', '_').replace('(', '').replace(')', '')
        rBAN_path = output + '/{}_peptideGraph.json'.format(ID)
        # Run rBAN with discoveryMode
        call('java -jar {}/../external/rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {}/{}_ -discoveryMode'.format(os.path.dirname(os.path.abspath(__file__)), ID, smiles, output, ID), shell=True)
        
    else:

        rBAN_path = rBAN

    return rBAN_path
# Check correctness of PeptideSeq
def check_pept_seq(ps):
    """
    The function checks the correctness of the processed structure. If we have 
    core peptide chain only from unrecognized monomers, we should break 
    BioCAT, because of it cant match nan structure on a clusters.
    In the case, when there are not any parsed monomers in the PeptideSeq excluding 'nan', 
    we return None and stop the precossing. In the case, when the PeptideSeq is 
    Variant1 = ['thr', 'pro'] and Variant2 = ['nan', 'nan']],
    Variant2 will be removed.

    Parameters
    ----------
    ps : dict
        Core peptide chains for different biosynthesis types (e.g. A, B, or C).
    Returns
    -------
    ps : dict
        Corrected core peptide chains for different biosynthesis pathways.
    None
        In the case, when the sequence contains only 'nan' monomers.
    """
    ps_cop = ps.copy()
    
    for bs_type in ps_cop:
        if ps_cop[bs_type] == {}:
            continue
            
        for var in ps[bs_type]:
            
            check = 0
            
            for tour in ps_cop[bs_type][var]:
                if tour  == ['nan'] * len(tour):
                    continue
                    
                check = 1

        if check == 0:
            
            ps[bs_type].pop(var)
    
    vals = list(ps.values())
    
    if vals == [{}] * len(vals): # If PeptideSeq is empty
        
        return None
    
    else:
        
        return ps

#Standartization of monomers
def make_standard(PeptideSeq):
    """
    The function change name of monomer, gives its by rBAN,
    on names from antiSMASH databases.

    Parameters
    ----------
    PeptideSeq : dict
        Core peptide chains for different biosynthesis types (e.g. A, B, or C).
    Returns
    -------
    PeptideSeq : dict
        Corrected core peptide chains for different biosynthesis types (e.g. A, B, or C).
    """
    hmms = os.path.dirname(os.path.abspath(__file__)) + '/../HMM/Bacteria_HMM'
    monomers_names_tax = listdir(hmms) #standardized list of possible to compare substrates 
    monomers_names = []

    for subst in monomers_names_tax:

        monomers_names.append(subst.split('_')[1][: -4])

    #For every variants getting standard monomer names
    for bios_path in PeptideSeq:
        if len(PeptideSeq[bios_path]) == 0:
            continue

        for var in PeptideSeq[bios_path]:
            for seq in PeptideSeq[bios_path][var]:
                for sub_AS in monomers_names:
                    for sub in seq:
                        if 'prob' in sub or 'prop' in sub:
                            continue

                        if 'dhpg' in sub:

                            PeptideSeq[bios_path][var][PeptideSeq[bios_path][var].index(seq)][seq.index(sub)] =  'dhpg' #Changing monomer name
                            continue
                        
                        if sub_AS in sub:
                            
                            PeptideSeq[bios_path][var][PeptideSeq[bios_path][var].index(seq)][seq.index(sub)] =  sub_AS #Changing monomer name

    for bios_path in PeptideSeq:
        if len(PeptideSeq[bios_path]) == 0:
            continue

        for var in PeptideSeq[bios_path]:
            for seq in PeptideSeq[bios_path][var]:
                for sub in seq:
                    if sub not in monomers_names:
                        
                        PeptideSeq[bios_path][var][PeptideSeq[bios_path][var].index(seq)][seq.index(sub)] =  'nan' #For substrates, which we have not HMM, giving nan name

    return PeptideSeq
# Make a dictionary with results to csv output
def get_results_to_csv(dict_res, output, id_smi):
    """
    Records Result file for processed substance.

    Parameters
    ----------
    dict_res : dict
        Results of BioCAT processing in dict format.
    output : str
        Path to BioCAT output directory.
    id_smi : str
        Name of the processed structure.
    """
    res_df = DataFrame(data=dict_res)
    res_df.to_csv('{}/Results_{}.tsv'.format(output, id_smi), sep='\t', index=False)           
