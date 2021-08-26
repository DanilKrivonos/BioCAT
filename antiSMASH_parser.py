import argparse 
from json import load
from pandas import DataFrame
from pandas import read_csv
from Get_AS_DF import get_df

#sortion genes in cluster and modules in genes
def cluster_sort(df):
    
    chromosomes = list(set(df['Name'].values))
    df2 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:

        Chr_df = df[df['Name'] == chromosome]
        genes = set(Chr_df['Gen ID'].values)
        #sorting genes in a BGC
        sort_indexes = []
        sort_indexes = list(df[df['Name'] == chromosome]['Start of gen'].sort_values().index)    
        sort_genes = []    
        sorted_df = df.iloc[sort_indexes]

        for ind in list(sorted_df['Gen ID'].values):
            if ind not in sort_genes:
                        
                sort_genes.append(ind)

        sort_indexes = {}
        #sorting modules in gene
        for gen in sort_genes:
            
            Gen_df = Chr_df[Chr_df['Gen ID'] == gen]
            sort_indexes[gen] = list(Gen_df['Protein start'].sort_values().index)

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

                    subcluster_id = '{}_{}'.format(cluster, sub)
                    df.loc[idx, 'ID'] = subcluster_id
    
    return df
#reversing negative directed genes
def reverse_neg(df):

    chromosomes = list(dict.fromkeys((df['Name'].values)))
    df1 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:

        Chr_df = df[df['Name'] == chromosome]
        clusters = list(dict.fromkeys(Chr_df['ID']))
        
        for bgc in clusters:
            
            BGC_df = Chr_df[Chr_df['ID'] == bgc]
            genes = list(dict.fromkeys(BGC_df['Gen ID']))
            SB = []
            genes_stack = [] 
            
            for gen in genes:

                Gen_df = BGC_df[BGC_df['Gen ID'] == gen]
                
                if '+' in list(Gen_df['Gen strand']):
                    if len(SB) != 0:
                        
                        genes_stack.extend(SB[-1: : -1])
                        SB = []
                        
                    genes_stack.append(gen)
                    
                else:
                    
                    SB.append(gen)
                    
            genes_stack.extend(SB[-1: : -1])
            
            for gen in genes_stack:

                Sort_gen_df = BGC_df[BGC_df['Gen ID'] == gen]
                df1 = df1.append(df.iloc[Sort_gen_df.index])

    df1.reset_index(drop=True, inplace=True)

    return df1
#************************THE FUNCTIONS RMOVING PART AFTER TE DOMAIN TO THE START OF CLUSTER************************
def shift_contig(df2, remove):
    for gen in remove:
        
        df2 = df2.append(gen)

    return df2
def remove_modificate(remove):
    
    remove2 = []
    deque = []
    first = 0
    first_gen = []
    
    for gen in remove:

        modules = set(gen.ModuleID.values)
        module_privi = ''
        
        for module in modules:

            Domain_BGC = gen[gen.ModuleID == module]
            domains = list(dict.fromkeys(Domain_BGC['Domain name'].values))
            
            for domain in domains:
                if 'Condensation' in domain or 'Cglyc' in domain:
                    
                    first = 1
            
        if first == 0:

            first_gen.append(gen)
            first = 0
            
        else:
            
            deque.append(gen)
            
    if len(first_gen):
        for gen in first_gen:
        
            remove2.append(gen)
        
    for gen in deque:
        
        remove2.append(gen)
        
    return remove2
                    
def sort_cluster_seq(df):
    
    chromosomes = list(dict.fromkeys(df['Name'].values))
    df2 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:

        Chr_df = df[df['Name'] == chromosome]
        clusters = list(dict.fromkeys(Chr_df['ID'].values))
        
        for bgc in clusters:

            remove = []
            stack_domain = []
            BGC_df = Chr_df[Chr_df['ID'] == bgc]
            genes = list(dict.fromkeys(BGC_df['Gen ID'].values))
            TE = 0
            gen_term = None
            gen_counter = 0 #Counting gen number to find index of terminating gene with TE domain
            
            for gen in genes:
                
                Gen_BGC = BGC_df[BGC_df['Gen ID'] == gen]
                modules = list(dict.fromkeys(Gen_BGC.ModuleID.values))
                module_privi = ''
                strand = list(set(Gen_BGC['Gen strand'].values))[0]
                
                for module in modules[:: -1]:

                    Domain_BGC = Gen_BGC[Gen_BGC.ModuleID == module]
                    domains = list(dict.fromkeys(Domain_BGC['Domain name'].values))
                    
                    if TE == 1:
      
                        remove.append(Domain_BGC)
                    
                    else:

                        stack_domain.append(Domain_BGC)
                    
                    for domain in domains:
                        if 'Thioesterase' in domain:

                            stack_domain.extend(remove)
                            remove = []
                            gen_term = gen_counter
                            TE = 1
                            
                gen_counter += 1
                
            if strand == '-' and not gen_term is None:
              
                last_gen = stack_domain[gen_term]
                stack_domain = stack_domain[: -1]
                remove = remove_modificate(remove)
                recombined = stack_domain 
                recombined += remove 
                recombined.append(last_gen)

            else:
                
                remove = remove_modificate(remove)
                recombined = remove + stack_domain
                #remove.extend(stack_domain)

            df2 = shift_contig(df2, recombined)
    return df2
#*******************************************************************************************************************************
def generate_table_from_antismash(json_path, table_out, dif_strand):
    
    out = table_out + '/table.tsv'
    get_df(json_path, out)
    df = read_csv(out, sep='\t')
    
    if dif_strand == None:

        subcluster = check_subcluster(df)

        if subcluster != None:

            df = split_subcluster(df, subcluster)

    df = cluster_sort(df)
    df = reverse_neg(df)
    df = sort_cluster_seq(df)
    df.to_csv(out, index=False, sep='\t')

print('Biosynthesis clusters have been successfully discovered!')
