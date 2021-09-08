import argparse 
from json import load
from pandas import DataFrame
from pandas import read_csv
from src.Get_AS_DF import get_df

""" This function construction correct NRPS clusters"""

# Sortion genes in cluster and modules in genes
def cluster_sort(df):
    """
    The function sorting clusters in table with meta inforamtion.
    Everything sorted from min to max coordinate, but it does not 
    correct in case of negative strand. Thus it will reversed by 
    reverse_neg function. Firstly sorting nucleotide sequence, 
    chromosomes and genes, then sortin domain by location in protein
    seques.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with NRPS meta information.
    Returns
    -------
    df2 : pandas DataFrame
        Corrected DataFrame with NRPS meta information.
    """
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
    """
    The function find genes with differing directions in area of one cluster.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with NRPS meta information.
    Returns
    -------
    subcluster : dict
        Dictionary  with information about potential subclusters.
    None
        If we have not potential subclusters.
    """
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
    """
    Split cluster on subclusters. In this case one cluster BGC_n can split
    on BGC_n_0, BGC_n_1, ...

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with NRPS meta information.
    subclusterdf : dict
        Dictionary  with information about potential subclusters.
    Returns
    -------
    df : pandas DataFrame
        Corrected DataFrame with NRPS meta information.
    """
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
    """
    Reversing negative strand genes from max to min coordinate.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with NRPS meta information.
    Returns
    -------
    df1 : pandas DataFrame
        Corrected DataFrame with NRPS meta information.
    """
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
    """
    The function append shifted fragment from 
    sort_cluster_seq function.

    Parameters
    ----------
    df2 : pandas DataFrame
        DataFrame NRPS cluster fragment.
    remove : list
        List of cluster fragment, which should removed.
    Returns
    -------
    df2 : pandas DataFrame
        Corrected DataFrame with NRPS meta information.
    """
    for gen in remove:
        
        df2 = df2.append(gen)

    return df2

def remove_modificate(remove):
    """
    If cluster has module with C starter domain it should
    be the first.

    Parameters
    ----------
    remove : list
        List of cluster fragment, which should removed.
    Returns
    -------
    dremove2 : list
        Corrected list of cluster fragment, which should removed.
    """
    remove2 = []
    deque = []
    first = 0
    first_gen = []
    
    for gen in remove:

        modules = set(gen.ModuleID.values)
        
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
    """
    In most cases the domain of Thioesterase sould be in the end of NRPS, 
    but it codding in the start of cluster. Thus, module with TE domain
    removed to the end.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame with NRPS meta information.
    Returns
    -------
    df2 : pandas DataFrame
        Corrected DataFrame with NRPS meta information.
    """
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
                remove = remove_modificate(remove) # Check first module
                recombined = stack_domain 
                recombined += remove 
                recombined.append(last_gen)

            else:
                
                remove = remove_modificate(remove)
                recombined = remove + stack_domain
                #remove.extend(stack_domain)
            # Recompiling of cluster
            df2 = shift_contig(df2, recombined)
    return df2
#*******************************************************************************************************************************
def generate_table_from_antismash(json_path, table_out, dont_dif_strand):
    """
    The function generate table with NRPS cluster mata inforamtion and 
    construct deque of modules of NRPS.

    Parameters
    ----------
    json_path : str
        Path to antiSMASH json output.
    table_out : str
        Path to BioCAT output directory.
    dont_dif_strand : bool
        Parameter of cutting cluster by directions of gene strands.
    """
    out = table_out + '/table.tsv'
    get_df(json_path, out)
    df = read_csv(out, sep='\t')
    
    if dont_dif_strand == False:

        subcluster = check_subcluster(df)

        if subcluster != None:

            df = split_subcluster(df, subcluster)

    df = cluster_sort(df)
    df = reverse_neg(df)
    df = sort_cluster_seq(df)
    df.to_csv(out, index=False, sep='\t')
