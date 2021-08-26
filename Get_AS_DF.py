from json import load 
from pandas import DataFrame

def get_info_dict(json_path):
    with open(json_path, 'r') as f:

        data = load(f)

    BIGDICT = {}

    for r in data['records']:

        dict_structure = {}
        dict_structure_proto = {}
        asModules = {}
        asmodule_id = 0
        target_CDSs = []
        asDomains = {}
        contig = r['id']
        types = []

        for feature in r['features']:

            loc = feature['location'].split

            if feature['type'] == 'cand_cluster':

                cand_id = feature['qualifiers']['candidate_cluster_number'][0]

                if cand_id not in dict_structure:
                    dict_structure[cand_id] = feature

            if feature['type'] == 'protocluster':
                #print(feature)
                proto_id = feature['qualifiers']['protocluster_number'][0]
                if proto_id not in dict_structure_proto:
                    dict_structure_proto[proto_id] = feature

            if feature['type'] == 'aSModule':
                target_CDSs.append(feature['qualifiers']['locus_tags'][0])
                asModules[asmodule_id] = feature
                asmodule_id += 1

            if feature['type'] == 'aSDomain':
                domain_id = feature['qualifiers']['domain_id'][0]
                asDomains[domain_id] = feature

            types.append(feature['type'])

        contig = r['id']

        for feature in r['features']:

            if feature['type'] != 'CDS':
                continue

            if feature['qualifiers']['locus_tag'][0] not in target_CDSs:
                continue

            start, end = feature['location'][1:-4].split(':')
            start, end = int(start), int(end)
            if 'locus_tag' in feature['qualifiers']:

                locus_tag = feature['qualifiers']['locus_tag'][0]

            elif 'protein_id' in feature['qualifiers']:

                locus_tag = feature['qualifiers']['protein_id'][0]

            if 'Module_childs' not in feature:

                feature['Module_childs'] = {}

            for asmodule in asModules:
                if asModules[asmodule]['qualifiers']['locus_tags'][0] == locus_tag:

                    feature['Module_childs'][asmodule] = asModules[asmodule]

                    for domain_id in asModules[asmodule]['qualifiers']['domains']:
                        if 'Domain_childs' not in feature['Module_childs'][asmodule]:

                            feature['Module_childs'][asmodule]['Domain_childs'] = {}
                        feature['Module_childs'][asmodule]['Domain_childs'][domain_id] = asDomains[domain_id]



            for protocluster in dict_structure_proto:

                loc = dict_structure_proto[protocluster]['location']
                pstart, pend = loc[1:-1].split(':')
                pstart, pend = int(pstart), int(pend)

                if pstart <= start <= end <= pend:
                    if 'CDS_childs' not in dict_structure_proto[protocluster]:

                        dict_structure_proto[protocluster]['CDS_childs'] = {}

                    dict_structure_proto[protocluster]['CDS_childs'][locus_tag] = feature

        for c in dict_structure:

            dict_structure[c]['childs'] = {}

            for protocluster in dict_structure[c]['qualifiers']['protoclusters']:

                dict_structure[c]['childs'][protocluster] = dict_structure_proto[protocluster]    

        BIGDICT[contig] = dict_structure
        
    return BIGDICT

def get_coord(gen_location):
    if 'join' in gen_location:
        print('nippnpnponopnopnop')
        if '-' in gen_location:

            start = gen_location.split(', ')[-1].split(':')[0][1: ]
            end = gen_location.split(', ')[0].split(':')[1].split(']')[0]
            strand = '-'
            
        elif '+' in gen_location:

            start = gen_location.split(', ')[0].split(':')[0].split('[')[1]
            end = gen_location.split(', ')[-1].split(':')[1].split(']')[0]
            strand = '+'
    else:

        start = gen_location.split(':')[0].split('[')[1]
        end = gen_location.split(':')[1].split(']')[0]
        strand = gen_location[-2]
        
    return start, end, strand

def get_df(json_path, out):
    
    BIGDICT = get_info_dict(json_path)
    cluster_id = 0
    protocluster_id = 0
    keys = {'Name' : [],
        'ID' : [],
        'Gen ID' : [],
        'Coordinates of cluster' : [],
        'Coordinates of protocluster' : [],
        'Gen strand' : [],
        'Start of gen' : [],
        'End of gen' : [],
        'Protein start' : [],
        'Protein end' : [],
        'ModuleID' : [],
        'Domain name' : [],
        'Sequence': []
        }

    for ID in BIGDICT.keys():

        Name = BIGDICT[ID]

        for cand_cluster in Name:

            BGC =  Name[cand_cluster]

            if 'NRPS' not in BGC['qualifiers']['product']:
                continue

            cluster_loc = BGC['location']
            
            
            for proto_cluster in BGC['childs']:

                proto_BGC = BGC['childs'][proto_cluster]
                proclust_loc = proto_BGC['location']
                
                if 'NRPS' not in proto_BGC['qualifiers']['product'][0]:
                    continue

                for gen in proto_BGC['CDS_childs'].keys():

                    gen_id = proto_BGC['CDS_childs'][gen]
                    gen_location = gen_id['location']
                    gen_start, gen_end, gen_strand = get_coord(gen_location)
                    mod_num = 1

                    for module in gen_id['Module_childs']:
                        
                        modlue_name = gen_id['Module_childs'][module]
                        module_loc = modlue_name['location']

                        for domain in modlue_name['Domain_childs']:

                            domain_info = modlue_name['Domain_childs'][domain]
                            location = domain_info['location']
                            protein_start = domain_info['qualifiers']['protein_start'][0]
                            protein_end = domain_info['qualifiers']['protein_end'][0]
                            translation = domain_info['qualifiers']['translation'][0]
                            # Protocluster
                            keys['Name'].append(ID)
                            keys['ID'].append('BGC_proto_{}_{}'.format(cluster_id, protocluster_id))
                            keys['Gen ID'].append(gen)
                            keys['Coordinates of cluster'].append(cluster_loc)
                            keys['Coordinates of protocluster'].append(proclust_loc)
                            keys['Gen strand'].append(gen_strand)
                            keys['Start of gen'].append(gen_start)
                            keys['End of gen'].append(gen_end)
                            keys['Protein start'].append(protein_start)
                            keys['Protein end'].append(protein_end)
                            keys['ModuleID'].append('{}_{}'.format(gen, mod_num))
                            keys['Domain name'].append(domain)
                            keys['Sequence'].append(translation)
                            # Candidate custer
                            keys['Name'].append(ID)
                            keys['ID'].append('BGC_cand_{}'.format(cluster_id))
                            keys['Gen ID'].append(gen)
                            keys['Coordinates of cluster'].append(cluster_loc)
                            keys['Coordinates of protocluster'].append(proclust_loc)
                            keys['Gen strand'].append(gen_strand)
                            keys['Start of gen'].append(gen_start)
                            keys['End of gen'].append(gen_end)
                            keys['Protein start'].append(protein_start)
                            keys['Protein end'].append(protein_end)
                            keys['ModuleID'].append('{}_{}'.format(gen, mod_num))
                            keys['Domain name'].append(domain)
                            keys['Sequence'].append(translation)
                    mod_num += 1
                protocluster_id += 1
            cluster_id += 1
            
    df = DataFrame(data=keys) 
    df.to_csv(out, index=False, sep='\t')