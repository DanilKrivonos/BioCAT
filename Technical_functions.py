import os
from json import load
from subprocess import call
from pandas import DataFrame
#Getting smiles formula and names of substances
def parse_smi_file(smiles, name, file_smiles):

    smile_list = []
    ids = []
    
    if smiles is not None:
        
        smile_list = [smiles]
        ids = [name.replace(' ', '_').replace('(', '').replace(')', '')]

    elif file_smiles is not None:
        with open(file_smiles) as smi:
            for line in smi:

                smile_list.append(line.split('\t')[1].replace('\n', ''))
                ids.append(line.split('\t')[0].replace(' ', '_').replace('(', '').replace(')', ''))
    
    return smile_list, ids

# Run antiSMASH if it neccecery or return path to directory 
def run_antiSMASH(antismash, output, genome, cpu):
    if antismash is None:
    
        anti_out = output + '/antismash_result/'

        try:
            
            os.mkdir(anti_out)

        except FileExistsError:

            print('The output directory already exists')
        # Run antiSMASH with loose mode
        call('antismash {} --cb-general --output-dir {} --genefinding-tool prodigal --cpus {}'.format(genome, anti_out, cpu), shell=True)
        json_path = anti_out + ('.').join(os.path.split(genome)[1].split('.')[0: -1]) + '.json'

    else:

        json_path = antismash

    return json_path

# Getting substance name if we have only rBAN json
def get_ids(outp):
    with open(outp, 'r') as json:

        js = load(json)
        
    try:

        name = js['id']

    except:

        name = js[0]['id']

    return  name
    
# Run rBAN if it neccecery or return path to directory 
def run_rBAN(rBAN, ID, SMI, output):
    #run rBUN
    print('Hydrolizing of substrate with rBAN ...')
    if rBAN is None:

        smiles = '"' + SMI + '"'
        #idsmi = ID.replace(' ', '_').replace('(', '').replace(')', '')
        rBAN_path = output + '/{}_peptideGraph.json'.format(ID)
        # Run rBAN with discoveryMode
        call('java -jar rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {}/{}_ -discoveryMode'.format(ID, smiles, output, ID), shell=True)
        
    else:

        rBAN_path = rBAN

    return rBAN_path
# Check correctness of PeptideSeq
def check_pept_seq(ps):
    
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
# Make a dictionary with results to csv output
def get_results_to_csv(dict_res, output, id_smi):

    res_df = DataFrame(data=dict_res)
    res_df.to_csv('{}/Results_{}.csv'.format(output, id_smi), sep='\t', index=False)           
