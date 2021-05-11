import argparse 
from json import load

def anti_parse(path, table_out):
    print(path)
    out = table_out + 'table.tsv'

    print('Starting parsing antiSMASH output ...')
    with open(out, 'w') as save_file:

        save_file.write('ID\tName\tCoordinates of cluster\tStrain\tCoordinates of module\tDomain name\tLarge prediction\tSmall prediction\tSingle aa prediction\tE-value\tScore\tSequence\n')

        with open(path, 'r') as js:

            translates = {}
            open_js = load(js)
            coord_tag = {}
            
            for key in open_js['records']:
                for i in key['features']:

                    if 'aSModule' not in i['type']:
                        continue

                    if 'nrps' not in i['qualifiers']['type']:
                        continue

                    location = i['location']
                    tag = i['qualifiers']['locus_tags'][0]
                    BGC = ''

                    for z in key['features']:

                        if 'aSDomain' not in z['type']:
                            continue

                        if z['qualifiers']['domain_id'][0] not in i['qualifiers']['domains']:
                            continue

                        evalue = z['qualifiers']['evalue'][0]
                        score = z['qualifiers']['score'][0]
                        translates[z['qualifiers']['domain_id'][0]] = z['qualifiers']['translation'] 

                    for dom in key['modules']['antismash.modules.nrps_pks']['domain_predictions']:
                        if dom in i['qualifiers']['domains']:
                            if 'NRPSPredictor2' in key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]:
                                
                                LargePred = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['large_cluster_pred']
                                SmalPred = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['small_cluster_pred']
                                Single_AA = key['modules']['antismash.modules.nrps_pks']['domain_predictions'][dom]['NRPSPredictor2']['single_amino_pred']
                            
                            else:
                                LargePred = 'nan'
                                SmalPred = 'nan'
                                Single_AA = 'nan'
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


                    if BGC == '':
                        continue
                    
                    for dom in key['modules']['antismash.modules.nrps_pks']['domain_predictions']:
                        if dom in i['qualifiers']['domains']:

                            save_file.write('{ID}\t{Name}\t{CL}\t{Str}\t{loc}\t{domain}\t{lp}\t{smap}\t{sip}\t{eva}\t{score}\t{translate}\n'.format(ID = BGC,
                                                                                                                                                    Name = key['id'].split('.')[0],
                                                                                                                                                    CL =  coord_tag,
                                                                                                                                                    Str = strand,
                                                                                                                                                    loc = location,
                                                                                                                                                    domain = dom,
                                                                                                                                                    lp = LargePred,
                                                                                                                                                    smap = SmalPred,
                                                                                                                                                    sip = Single_AA,
                                                                                                                                                    eva = evalue,
                                                                                                                                                    score = score,
                                                                                                                                                    translate = translates[dom][0]))
print('Biosynthesis clusters have been successfully discovered!')