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
#************************THIS TWO FUNCTION RMOVING PART AFTER TE DOMAIN TO THE START OF CLUSTER************************
def shift_contig(df2, df, remove):
    for domain in remove:

        Dom_df = df[df['Domain name'] == domain]
        df2 = df2.append(Dom_df)
        
    return df2

def sort_cluster_seq(df):
    
    chromosomes = list(dict.fromkeys(df['Name'].values))
    df2 = DataFrame(columns=list(df.keys()))

    for chromosome in chromosomes:

        Chr_df = df[df['Name'] == chromosome]
        clusters = list(dict.fromkeys(Chr_df['ID']))
        
        for bgc in clusters:
            
            add = 0
            remove = []
            stack_domain = []
            BGC_df = Chr_df[Chr_df['ID'] == bgc]
            domains = BGC_df['Domain name']
            
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
