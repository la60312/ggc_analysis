import pandas as pd
import numpy as np
import os
import json
import pprint
import argparse

def is_variant_found(data):
    is_found = False
    results = data['test_result']
    if 'result' in results:
        r = results['result']
        if r == 'VARIANTS_DETECTED':
            is_found = True
    return is_found

def is_exome_found(data):
    is_found = False
    if 'test_type' in data:
        results = data['test_type']
        is_found = results == 'EXOME_SEQUENCING'
    return is_found

def is_pathogenic(v):
    found = False
    var_class = ''
    if 'mutation' in v and v['mutation'] != []:
        var_class = v['mutation']['interpretation']
        if var_class == 'PATHOGENIC' or var_class == 'LIKELY_PATHOGENIC' or var_class == 'UNCERTAIN_SIGNIFICANCE':
            found = True
    if 'mutation1' in v:
        var_class = v['mutation1']['interpretation']
        if var_class == 'PATHOGENIC' or var_class == 'LIKELY_PATHOGENIC' or var_class == 'UNCERTAIN_SIGNIFICANCE':
            found = True
    return found, var_class


parser = argparse.ArgumentParser(description='Get summary table of results')
parser.add_argument('-i', '--input', help='Path of input')
parser.add_argument('-o', '--output', default=".", help='Path of output')
parser.add_argument('-g', '--genomic', default=".", help='Path of genomic data')
args = parser.parse_args()

result_dir = args.output

p_dir = args.input

# The genomic info of GGC cases
genomic_file = args.genomic
with open(genomic_file, 'rb') as f:
    genomic_data = json.load(f)

genomics = {}
for g in genomic_data:
    if 'case_id' in g['genomic_data']:
        case_id = str(g['genomic_data']['case_id'])
        if case_id in genomics:
            genomics[case_id].append(g)
        else:
            genomics[case_id] = [g]

# types of using different combination of scores
test_types = ['pedia', 'f_g', 'f_c_g', 'c_p_b', 'c', 'g', 'f']

# output dataframe
out_df = pd.DataFrame()
flag = True
for t_idx, test_type in enumerate(test_types):
    print('----------------------------------------')
    # Init
    input_dir = os.path.join(p_dir, 'ggc_' + test_type)
    input_results_files = [f for f in os.listdir(input_dir) if 'json' in f]
    case_list = []
    all_match_gene = []
    all_match_rank = []
    all_interp = []
    all_c = []
    all_g = []
    all_f = []
    all_p = []
    all_b = []
    all_syn = []
    all_diag = []
    all_mask = []
    all_exome_found = []
    all_no_found = []
    all_features = []

    # Iterate all file in input folder
    for input_results_file in input_results_files:
        # Parse patient JSON file
        case_id = input_results_file[0:-5]
        with open(os.path.join(input_dir, input_results_file), 'r') as f:
            j_data = json.load(f)
            syn_data = j_data['selected_syndromes']
            syns = [syn['syndrome_name'] for syn in syn_data]
            diag = [syn['diagnosis'] for syn in syn_data]
            mask = [syn['has_mask'] for syn in syn_data]

        # Read result file
        with open(os.path.join(input_dir, input_results_file[0:-5]+'.csv'), 'r') as f:
            result = pd.read_csv(f)

        # PEDIA cases don't have the same genomic entry format
        pedia_case = False
        if case_id not in genomics:
            if len(j_data['genomicData']) == 0:
                #all_no_found.append(False)
                print('case {} genomic data not found'.format(case_id))
                continue
            # if genomicData is not empty, this case is from PEDIA real case
            exome_found = [True]
            pedia_case = True
            all_no_found.append(True)
            found_data = []
        else:
            g_datas = genomics[case_id]
            exome_found = [is_exome_found(g['genomic_data']) for g in g_datas]
            found_data = [g for g in g_datas if is_variant_found(g['genomic_data'])]
            all_no_found.append(True)

        # Append meta data
        all_features.append(j_data['features'])
        all_syn.append(str(syns)[1:-1])
        all_diag.append(str(diag)[1:-1])
        all_mask.append(str(mask)[1:-1])

        # Check if exome test is found
        x = True in exome_found
        all_exome_found.append(x)

        match_gene = []
        match_rank = []
        interp = []
        c = []
        ge = []
        f = []
        p = []
        b = []

        # By default, we assume the case is from GGC, so we parse
        # the genomic json file to get gene info
        # and further get rank and scores
        for g in found_data:
            if case_id == '62094':
                print(g)
            variants = g['genomic_data']['test_result']['variants']
            for v in variants:
                #print(v)
                if 'gene' not in v:
                    continue
                is_p, v_class = is_pathogenic(v)
                if not is_p:
                    continue

                # Sometimes gene_symbol is empty, have to check gene_name as well
                gene_id = v['gene']['gene_symbol']
                if gene_id == "":
                    gene_id = v['gene']['gene_name']

                # Find the rank of correct disease-causing gene
                index = result.index[result['gene_name'] == gene_id]

                # If gene is not found in PEDIA results
                if len(index) == 0:
                    print('{} variant not found'.format(case_id))
                    if gene_id not in match_gene:
                        match_gene.append(gene_id)
                        interp.append(v_class)
                        match_rank.append('Not found')
                else:
                    if gene_id not in match_gene:
                        match_gene.append(gene_id)
                        match_rank.append(index[0] + 1)
                        ge.append(result['gestalt_score'][index[0]])
                        f.append(result['feature_score'][index[0]])
                        c.append(result['cadd_score'][index[0]])
                        p.append(result['pheno_score'][index[0]])
                        b.append(result['boqa_score'][index[0]])
                        interp.append(v_class)

        # If it is PEDIA cases, then we parse the genomic data in patient JSON directly
        if pedia_case:
            # the disease-causing gene
            gene_id = j_data['genomicData'][0]['Test Information']['Gene Name']

            index = result.index[result['gene_name'] == gene_id]
            match_gene.append(gene_id)
            match_rank.append(index[0] + 1)
            interp.append('PATHOGENIC')
            ge.append(result['gestalt_score'][index[0]])
            f.append(result['feature_score'][index[0]])
            c.append(result['cadd_score'][index[0]])
            p.append(result['pheno_score'][index[0]])
            b.append(result['boqa_score'][index[0]])

        all_match_rank.append(str(match_rank)[1:-1])
        all_match_gene.append(str(match_gene)[1:-1])
        all_interp.append(str(interp)[1:-1])
        all_g.append(str(ge)[1:-1])
        all_f.append(str(f)[1:-1])
        all_c.append(str(c)[1:-1])
        all_p.append(str(p)[1:-1])
        all_b.append(str(b)[1:-1])
        case_list.append(case_id)

    if flag:
        out_df['case_id'] = case_list
        flag = False
        out_df['exome test found'] = all_exome_found
        out_df['test found'] = all_no_found
        out_df['syndrome'] = all_syn
        out_df['num of features'] = [len(fs) for fs in all_features]

        out_df['diagnosis'] = all_diag
        out_df['mask'] = all_mask
        out_df['gene'] = all_match_gene
        out_df['classification'] = all_interp
        out_df['gestalt_score'] = all_g
        out_df['fm_score'] = all_f
        out_df['cadd_score'] = all_c
        out_df['pheno_score'] = all_p
        out_df['boqa_score'] = all_b

    out_df[test_type+'_rank'] = all_match_rank
    if t_idx == len(test_types)-1:
        out_df['features'] = all_features

out_df.to_csv(os.path.join(args.output, 'ggc_results.csv'), header=True, index=False)
