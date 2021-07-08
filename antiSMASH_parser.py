import argparse 
from json import load
from pandas import DataFrame

#sortion genes in cluster and modules in genes
def cluster_sort(df):
    
    chromosomes = list(set(df['Name'].values))
    df2 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:

        genes = set(df[df['Name'] == chromosome]['Gen ID'].values)
        #sorting genes in a BGC
        sort_indexes = []
        sort_indexes = list(df[df['Name'] == chromosome]['Start of gen'].sort_values().index)    
        sort_genes = []    
        
        for ind in list(df.iloc[sort_indexes]['Gen ID'].values):
            if ind not in sort_genes:
                        
                sort_genes.append(ind)

        sort_indexes = {}
        #sorting modules in gene
        for gen in sort_genes:
                
            sort_indexes[gen] = list(df[df['Gen ID'] == gen]['Protein start'].sort_values().index)

        for idx in sort_indexes:

            df2 = df2.append(df.iloc[sort_indexes[idx]])
    df2.reset_index(drop=True, inplace=True)
    return df2
#finding subclusters
def check_subcluster(df):

    chromosomes = list(set(df['Name'].values))
    have_not = 0
    subcluster = {}
    
    for chromosome in chromosomes:
        
        cluster_stack = df[df['Name'] == chromosome]['ID'].values
        subcluster[chromosome] = {}

        for cluster in cluster_stack:
            
            subcluster[chromosome][cluster] = {}
            strand = df[df['ID'] == cluster]['Gen strand'].values[0]
            cluster_stack = set(cluster_stack)

            if all(gen == strand for gen in df[df['ID'] == cluster]['Gen strand'].values):
                
                continue
                
            else:
                
                have_not = 1
                subcluster[chromosome][cluster] = {0: [],
                                                   1: []} 

                plus_SB = []
                minus_SB = []
        
                for idx in df[df['ID'] == cluster].index:
                    if '+' in df['Gen strand'][idx]:
                        
                        plus_SB.append(idx)
                        
                    else:
                        
                        minus_SB.append(idx)
                        
                subcluster[chromosome][cluster][0].extend(plus_SB)
                subcluster[chromosome][cluster][1].extend(minus_SB)

    if have_not != 0:

        return subcluster
    
    else:
        
        return None #if cluster have not diffetent direction
#Spliiting subclusters 
def split_subcluster(df, subcluster):

    chromosomes = list(set(df['Name'].values))

    for chromosome in chromosomes:
        for cluster in subcluster[chromosome].keys():
            if len(subcluster[chromosome][cluster]) < 2:
                continue
                
            subclusters = subcluster[chromosome][cluster]

            for sub in subclusters.keys():
            
                for idx in df.iloc[subclusters[sub]].index:
                    
                    df['ID'][idx] = '{}_{}'.format(cluster, sub)
    
    return df
#reversing negative directed genes
def reverse_neg(df):

    chromosomes = list(dict.fromkeys((df['Name'].values)))
    df1 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:
        
        clusters = list(dict.fromkeys(df[df['Name'] == chromosome]['ID']))
        
        for bgc in clusters:
            
            genes = list(dict.fromkeys(df[df['ID'] == bgc]['Gen ID']))
            SB = []
            genes_stack = [] 
            
            for gen in genes:
                
                if '+' in list(df[df['ID'] == bgc][df['Gen ID'] == gen]['Gen strand']):
                    if len(SB) != 0:
                        
                        genes_stack.extend(SB[-1: : -1])
                        SB = []
                        
                    genes_stack.append(gen)
                    
                else:
                    
                    SB.append(gen)
                    
            genes_stack.extend(SB[-1: : -1])
            
            for gen in genes_stack:

                df1 = df1.append(df.iloc[df[df['ID'] == bgc][df['Gen ID'] == gen].index])

    df1.reset_index(drop=True, inplace=True)

    return df1
#************************THIS TWO FUNCTION RMOVING PART AFTER TE DOMAIN TO THE START OF CLUSTER************************
def shift_contig(df2, df, remove):
    for domain in remove:

        df2 = df2.append(df[df['Domain name'] == domain])
        
    return df2

def sort_cluster_seq(df):
    
    chromosomes = list(dict.fromkeys(df['Name'].values))
    df2 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:
        
        clusters = list(dict.fromkeys(df[df['Name'] == chromosome]['ID']))
        
        for bgc in clusters:
            
            add = 0
            remove = []
            stack_domain = []
            domains = df[df['ID'] == bgc]['Domain name']
            
            for domain in domains:
                if add == 1:
                    
                    remove.append(domain)
                    
                else:
                    
                    stack_domain.append(domain)
                    
                if 'Thioesterase' in domain:
                    
                    stack_domain.extend(remove)
                    remove = []
                    add = 1

            remove.extend(stack_domain)
            df2 = shift_contig(df2, df, remove)
    
    return df2
#*******************************************************************************************************************************
def anti_parse(path, table_out, dif_strand):
   
    out = table_out + '/table.tsv'

    print('Starting parsing antiSMASH output ...')
    keys = {'Name' : [],
            'ID' : [], 
            'Gen ID' : [],
            'Coordinates of cluster' : [],
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
        check_gen_id = {}
        chromosome_ID = 0

        for key in open_js['records']:
            for i in key['features']:
                if 'protocluster' in i['type']:
                    
                    BGC = 'BGC_{}_{}'.format(chromosome_ID, i['qualifiers']['protocluster_number'][0])
                    coord_cluster = i['location']
                
                if 'CDS' == i['type']:
                    try:
                        if 'locus_tag' in i['qualifiers']:

                            check_gen_id[i['qualifiers']['locus_tag'][0]] = i['qualifiers']['locus_tag'][0]
                            coordhave = 'locus_tag'
                            GEN_COORD[i['qualifiers']['locus_tag'][0]] = i['location']
                            
                        elif 'gene' in i['qualifiers']:
                            
                            gene_loc_tag[i['qualifiers']['gene'][0]] = i['qualifiers']['protein_id'][0]
                            coordhave = 'gene'
                            GEN_COORD[i['qualifiers']['protein_id'][0]] = i['location']
                            
                        elif 'protein_id' in i['qualifiers']:

                            coordhave = 'protein_id'
                            GEN_COORD[i['qualifiers']['protein_id'][0]] = i['location']
                            
                    except:
                        continue
                
                if 'aSModule' not in i['type']:
                    continue

                tag = i['qualifiers']['locus_tags'][0]
                metric_dict = {}

                for z in key['features']:

                    if 'aSDomain' not in z['type']:
                        continue

                    metric_dict[z['qualifiers']['domain_id'][0]] = {'evalue' : z['qualifiers']['evalue'][0],
                                                                    'score' : z['qualifiers']['score'][0]}
                    domain_name = z['qualifiers']['domain_id'][0]
                    protein_start[domain_name] = z['qualifiers']['protein_start'][0]
                    protein_end[domain_name] = z['qualifiers']['protein_end'][0]
                    translates[domain_name] = z['qualifiers']['translation']                         

                    if coordhave == 'gene':

                        GEN_ID[domain_name] = gene_loc_tag[z['qualifiers']['locus_tag'][0]]
                    
                    elif coordhave == 'locus_tag':
                        if z['qualifiers']['locus_tag'][0] not in check_gen_id:
                            continue
                            
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
                                        
                                predictiton_dict[dom]['Single_AA'] = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['stachelhaus_predictions'][0].replace('(', '').replace(')', '')

                        else:
                            
                            predictiton_dict[dom]['LargePred'] = 'nan'
                            predictiton_dict[dom]['SmalPred'] = 'nan'
                            predictiton_dict[dom]['Single_AA'] = 'nan'

                for dom in key['modules']['antismash.modules.nrps_pks']['domain_predictions']:
                    if dom in i['qualifiers']['domains']:
                        if dom not in GEN_ID:
                            continue
                            
                        if 'join' in GEN_COORD[GEN_ID[dom]]:
                            if '-' in GEN_COORD[GEN_ID[dom]]:

                                start = GEN_COORD[GEN_ID[dom]].split(', ')[-1].split(':')[0][1: ]
                                end = GEN_COORD[GEN_ID[dom]].split(', ')[0].split(':')[1].split(']')[0]

                            elif '+' in GEN_COORD[GEN_ID[dom]]:

                                start = GEN_COORD[GEN_ID[dom]].split(', ')[0].split(':')[0].split('[')[1]
                                end = GEN_COORD[GEN_ID[dom]].split(', ')[-1].split(':')[1].split(']')[0]
                            
                        else:
                            
                            start = GEN_COORD[GEN_ID[dom]].split(':')[0].split('[')[1]
                            end = GEN_COORD[GEN_ID[dom]].split(':')[1].split(']')[0]
                        
                        keys['Name'].append(key['id'].split('.')[0])    
                        keys['ID'].append(BGC)
                        keys['Gen ID'].append(str(GEN_ID[dom]))
                        keys['Coordinates of cluster'].append(coord_cluster)
                        keys['Gen strand'].append(str(GEN_strand[dom]))
                        keys['Start of gen'].append(int(start))
                        keys['End of gen'].append(int(end))
                        keys['Protein start'].append(int(protein_start[dom]))
                        keys['Protein end'].append(int(protein_end[dom]))
                        keys['Domain name'].append(dom)
                        keys['Large prediction'].append(predictiton_dict[dom]['LargePred'])
                        keys['Small prediction'].append(predictiton_dict[dom]['SmalPred'])
                        keys['Single aa prediction'].append(predictiton_dict[dom]['Single_AA'])
                        keys['E-value'].append(metric_dict[dom]['evalue'])
                        keys['Score'].append(metric_dict[dom]['score'])
                        keys['Sequence'].append(translates[dom][0])

            chromosome_ID += 1

    df = DataFrame(data=keys)

    if dif_strand == None:

        subcluster = check_subcluster(df)

        if subcluster != None:

            df = split_subcluster(df, subcluster)

    df.to_csv(out, index=False, sep='\t')
    from pandas import read_csv 
    df = read_csv(out, sep='\t')
    df = cluster_sort(df)
    df = reverse_neg(df)
    df = sort_cluster_seq(df)
    df.to_csv(out, index=False, sep='\t')

print('Biosynthesis clusters have been successfully discovered!')
