import os
from json import load
from subprocess import call

#Getting smiles formula and names of substances
def parse_smi_file(smiles, name, file_smiles):

    smile_list = []
    ids = []
    
    if smiles is not None:
        
        smile_list = [smiles]
        ids = [name]

    elif file_smiles is not None:
        with open(file_smiles) as smi:
            for line in smi:

                smile_list.append(line.split('\t')[1].replace('\n', ''))
                ids.append(line.split('\t')[0])
    
    return smile_list, ids

# Run antiSMASH if it neccecery or return path to directory 
def run_antiSMASH(antismash, output, genome, cpu):
    if antismash is None:
    
        anti_out = output + 'antismash_result/'

        try:
            
            os.mkdir(anti_out)

        except FileExistsError:

            print('The output directory already exists')
        # Run antiSMASH with loose mode
        call('antismash {} --cb-general --output-dir {} --genefinding-tool prodigal --hmmdetection-strictness loose --cpus {}'.format(genome, anti_out, cpu), shell=True)
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
        idsmi = ID.replace(' ', '_').replace('(', '').replace(')', '')
        rBAN_path = output + '/{}_peptideGraph.json'.format(idsmi)
        # Run rBAN with discoveryMode
        call('java -jar rBAN-1.0.jar -inputId {} -inputSmiles {} -outputFolder {}/{}_ -discoveryMode'.format(idsmi, smiles, output, idsmi), shell=True)
        
    else:

        rBAN_path = rBAN

    return rBAN_path