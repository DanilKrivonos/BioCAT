import argparse 
from json import load
from pandas import DataFrame

#sortion genes in cluster and modules in genes
def cluster_sort(df):
    
    genes = set(df['Gen ID'].values)
    #sorting genes in a BGC
    sort_indexes = []
    sort_indexes = list(df['Start of gen'].sort_values().index)    
    sort_genes = []    
    
    for ind in list(df.iloc[sort_indexes]['Gen ID'].values):
        if ind not in sort_genes:
                    
            sort_genes.append(ind)

    sort_indexes = {}
    df2 = DataFrame(columns=list(df.keys()))
    #sorting modules in gene
    for gen in sort_genes:
        
        sort_indexes[gen] = list(df[df['Gen ID'].str.contains(gen)]['Protein start'].sort_values().index)

    for idx in sort_indexes:
        
        df2 = df2.append(df.iloc[sort_indexes[idx]])

    return df2
#finding subclusters
def check_subcluster(df):
    
    cluster_stack = df['ID'].values
    cluster_stack = set(cluster_stack)
    subcluster = {}
    have_not = 0
    
    for cluster in cluster_stack:
        
        subcluster[cluster] = {}
        strand = df[df['ID'] == cluster]['Gen strand'].values[0]
        
        if all(gen == strand for gen in df[df['ID'] == cluster]['Gen strand'].values):
            
            continue
            
        else:
            
            genes = df[df['ID'] == cluster]['Gen strand']
            index = 0

            for gen_id in genes.index:

                if index not in subcluster[cluster]:
                    
                    subcluster[cluster][index] = []

                if genes[gen_id] != strand:
                    
                    index += 1
                    strand = genes[gen_id]
                    have_not = 1
                    
                    if index not in subcluster[cluster]:
                        
                        subcluster[cluster][index] = []
                    
                subcluster[cluster][index].append(gen_id)
    if have_not != 0:
        
        return subcluster
    
    else:
        
        return None #if cluster have not diffetent direction
#reversing negative directed genes
def reverse_neg(df):
    
    genes = list(dict.fromkeys(df['Gen ID'].values))
    gen_stack = {}
    df1 = DataFrame(columns=list(df.keys()))
    
    for gen in genes:
        if '-' in df[df['Gen ID'].str.contains(gen)]['Gen strand'].values:
            print(df[df['Gen ID'].str.contains(gen)]['Domain name'])
            print(list(df[df['Gen ID'].str.contains(gen)].index))
            gen_stack[gen] = list(df[df['Gen ID'].str.contains(gen)].index)
 

        if '+' in df[df['Gen ID'].str.contains(gen)]['Gen strand'].values:
            if gen_stack != {}:
                for gen_k in list(gen_stack.keys())[: : -1]:
                   
                    dfcop = df.iloc[gen_stack[gen_k]]
                    df1 = df1.append(dfcop)                
            
            dfcop = df[df['Gen ID'].str.contains(gen)]
            df1 = df1.append(dfcop)
            gen_stack = {}              
            
    if gen_stack != {}:
        for gen_k in list(gen_stack.keys())[: : -1]:

            dfcop = df.iloc[gen_stack[gen_k]]
            df1 = df1.append(dfcop)   
            
    return df1
#Spliiting subclusters 
def split_subcluster(df, subcluster):
    for cluster in subcluster.keys():
        if len(subcluster[cluster]) < 2:
            continue
            
        subclusters = subcluster[cluster]
        
        for sub in subclusters.keys():
        
            for idx in df.iloc[subclusters[sub]].index:
                df['ID'][idx] = '{}_{}'.format(cluster, sub)
    
    return df

def anti_parse(path, table_out):

    out = table_out + 'table.tsv'
    print('Starting parsing antiSMASH output ...')
    keys = {'ID' : [], 
            'Gen ID' : [],
            'Name' : [],
            'Coordinates of cluster' : [],
            'Strand' : [],
            'Gen strand' : [],
            'Start of gen' : [],
            'End of gen' : [],
            'Protein start' : [],
            'Protein end' : [],
            'Domain name' : [],
            'Large prediction' : [],
            'Small prediction' : [],
            'Single aa prediction' : [],
            'E-value' : [],
            'Score' : [],
            'Sequence': []
            }

    with open(path, 'r') as js:

        translates = {}
        open_js = load(js)
        coord_tag = {}
        protein_start = {}
        protein_end = {}
        GEN_ID = {}
        GEN_strand = {}
        GEN_COORD = {}
        gene_loc_tag = {}
        coordhave = ''
        chack_gen_id = {}

        for key in open_js['records']:
            for i in key['features']:
                if 'CDS' == i['type']:
                    try:
                        if 'gene' in i['qualifiers']:
                            
                            coordhave = 'gene'
                            gene_loc_tag[i['qualifiers']['gene'][0]] = i['qualifiers']['protein_id'][0]

                        elif 'locus_tag' in i['qualifiers']:

                            check_gen_id[i['qualifiers']['locus_tag'][0]] = i['qualifiers']['locus_tag'][0]

                            if i['qualifiers']['locus_tag'][0] not in GEN_COORD:
                                
                                coordhave = 'locus_tag' 
                                GEN_COORD[i['qualifiers']['locus_tag'][0]] = {'start' : i['location'].split('(')[0].split(':')[0][1: ],
                                                                             'end' : i['location'].split('(')[0].split(':')[1][: -1]}

                            if 'join' not in i['location']:

                                GEN_COORD[i['qualifiers']['locus_tag'][0]]['start'] = i['location'].split('(')[0].split(':')[0][1: ]
                                GEN_COORD[i['qualifiers']['locus_tag'][0]]['end'] = i['location'].split('(')[0].split(':')[1][: -1]

                            else:

                                if '-' in i['location']:

                                    start = i['location'].split(', ')[-1].split(':')[0][1: ]
                                    end = i['location'].split(', ')[0].split(':')[1].split(']')[0]
                                    GEN_COORD[i['qualifiers']['locus_tag'][0]]['start'] = start
                                    GEN_COORD[i['qualifiers']['locus_tag'][0]]['end'] = end

                                elif '+' in i['location']:

                                    start = i['location'].split(', ')[0].split(':')[0].split('[')[1]
                                    end = i['location'].split(', ')[-1].split(':')[1].split(']')[0]
                                    GEN_COORD[i['qualifiers']['locus_tag'][0]]['start'] = start
                                    GEN_COORD[i['qualifiers']['locus_tag'][0]]['end'] = end
                                    
                        elif 'protein_id' in i['qualifiers']:
                            
                            coordhave = 'protein_id'
                            GEN_COORD[i['qualifiers']['protein_id'][0]] = {'start' : i['location'].split(':')[0][1: ],
                                                                             'end' : i['location'].split(':')[1][: -1]}
                            
                    except:
                        continue

                if 'aSModule' not in i['type']:
                    continue

                if 'nrps' not in i['qualifiers']['type']:
                    continue

                tag = i['qualifiers']['locus_tags'][0]
                BGC = ''
                metric_dict = {}
                
                for z in key['features']:

                    if 'aSDomain' not in z['type']:
                        continue

                    if z['qualifiers']['domain_id'][0] not in i['qualifiers']['domains']:
                        continue

                    metric_dict[z['qualifiers']['domain_id'][0]] = {'evalue' : z['qualifiers']['evalue'][0],
                                                                    'score' : z['qualifiers']['score'][0]}
                    domain_name = z['qualifiers']['domain_id'][0]
                    protein_start[domain_name] = z['qualifiers']['protein_end'][0]
                    protein_end[domain_name] = z['qualifiers']['protein_start'][0]
                    translates[domain_name] = z['qualifiers']['translation'] 
                    
                    if coordhave == 'gene':
                        
                        GEN_ID[domain_name] = gene_loc_tag[z['qualifiers']['locus_tag'][0]]
                        
                        if z['qualifiers']['locus_tag'][0] not in GEN_COORD:

                            GEN_COORD[gene_loc_tag[z['qualifiers']['locus_tag'][0]]] = {'start' :'',
                                                                                     'end' : ''}
                    elif coordhave == 'loqus_tag':

                        GEN_ID[domain_name] = check_gen_id[z['qualifiers']['locus_tag'][0]]
                        
                    elif coordhave == 'protein_id':
                        
                        if z['qualifiers']['locus_tag'][0] in GEN_COORD:
                            
                            GEN_ID[domain_name] = z['qualifiers']['locus_tag'][0]
                        
                    if '+' in z['location']:
                        
                        direct = '+'
                        
                    elif '-' in z['location']:
                        
                        direct = '-'
                        
                    GEN_strand[domain_name] = direct
                predictiton_dict = {}
                
                for dom in key['modules']['antismash.modules.nrps_pks']['domain_predictions']:
                    
                    predictiton_dict[dom] = {'LargePred' : '',
                                            'SmalPred' : '',
                                            'Single_AA' : ''}
                    
                    if dom in i['qualifiers']['domains']:
                        if 'NRPSPredictor2' in key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]:
                            predictiton_dict[dom]['LargePred'] = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['large_cluster_pred']
                            predictiton_dict[dom]['SmalPred'] = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['small_cluster_pred']
                            predictiton_dict[dom]['Single_AA'] = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['single_amino_pred']
                            
                            if predictiton_dict[dom]['Single_AA'] == 'N/A':
                                        
                                        predictiton_dict[dom]['Single_AA'] = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['stachelhaus_predictions'][0]
                        else:
                            
                            predictiton_dict[dom]['LargePred'] = 'nan'
                            predictiton_dict[dom]['SmalPred'] = 'nan'
                            predictiton_dict[dom]['Single_AA'] = 'nan'
                            
                cutoff = 0

                for knw in  key['modules']['antismash.detection.hmm_detection']['rule_results']['cds_by_protocluster']:
                    for knw_i in knw[1]:
                        if tag in knw_i['cds_name']:

                            if cutoff < eval(knw[0]['qualifiers']['cutoff'][0]):

                                cutoff = eval(knw[0]['qualifiers']['cutoff'][0])
                                coord_tag = knw[0]['location']

                for knw in  key['modules']['antismash.modules.clusterblast']['knowncluster']['results']:
                    for suknw in knw['ranking']:
                        if 'NRP' not in suknw[0]['cluster_type']:
                            continue 
                            
                        if BGC != '':
                            break
                            
                        for i_k in suknw[1]['pairings']:

                            if tag == i_k[0].split('|')[-2]:
                                
                                BGC = suknw[0]['accession']
                                strand = i_k[-1]['strand']
                                
                if coordhave == 'gene' or coordhave == 'protein_id':
                    for knw in  key['modules']['antismash.modules.clusterblast']['knowncluster']['proteins']:
                        if knw['name'] in GEN_COORD:

                            GEN_COORD[knw['name']]['start'] = knw['location'].split('-')[0]
                            GEN_COORD[knw['name']]['end'] = knw['location'].split('-')[1]

                if BGC == '':
                    continue

                for dom in key['modules']['antismash.modules.nrps_pks']['domain_predictions']:
                    if dom in i['qualifiers']['domains']:
                        
                        keys['ID'].append(BGC)
                        keys['Gen ID'].append(GEN_ID[dom])
                        keys['Name'].append(key['id'].split('.')[0])
                        keys['Coordinates of cluster'].append(coord_tag)
                        keys['Strand'].append(strand)
                        keys['Gen strand'].append(GEN_strand[dom])
                        keys['Start of gen'].append(int(GEN_COORD[GEN_ID[dom]]['start']))
                        keys['End of gen'].append(int(GEN_COORD[GEN_ID[dom]]['end']))
                        keys['Protein start'].append(int(protein_start[dom]))
                        keys['Protein end'].append(int(protein_end[dom]))
                        keys['Domain name'].append(dom)
                        keys['Large prediction'].append(predictiton_dict[dom]['LargePred'])
                        keys['Small prediction'].append(predictiton_dict[dom]['SmalPred'])
                        keys['Single aa prediction'].append(predictiton_dict[dom]['Single_AA'])
                        keys['E-value'].append(metric_dict[dom]['evalue'])
                        keys['Score'].append(metric_dict[dom]['score'])
                        keys['Sequence'].append(translates[dom][0])

    df = DataFrame(data=keys)
    df = cluster_sort(df)
    subcluster = check_subcluster(df)

    if subcluster != None:

        df = split_subcluster(df, subcluster)
    
    df.to_csv(out, index=False, sep='\t')
    from pandas import read_csv 
    df = read_csv(out, sep='\t')
    df = reverse_neg(df)
    df.to_csv(out, index=False, sep='\t')
print('Biosynthesis clusters have been successfully discovered!')
