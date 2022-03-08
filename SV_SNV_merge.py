import pandas as pd
import sys
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-dfolder", help="put all sv and snv data you want to merge in the folder, same sample merging files must have same para and in  *SV* *SNV* pattern")
parser.add_argument("-g", help="gene list in csv formate")
parser.add_argument("-o", help="output file name")
args = parser.parse_args()

# data import
sv_data_list = glob.glob(f'{args.dfolder}/*SV*')
snv_data_list = glob.glob(f'{args.dfolder}/*SNV*')
sv_data_list.sort()
snv_data_list.sort()
gene = pd.read_csv(args.g, sep='\t', header=None)



def sv_data_arrange(sv_data, num):
    data = pd.read_csv(sv_data, sep=',')
    col = ['Chr', 'Start', 'End', 'SV_length', 'SV_type', 'Gene_name', 'Func.refGene/Location1','ExonicFunc.refGene/location2', \
           'Variant_ID', 'GD_AF_EAS', 'GD_POPMAX_AF', 'ACMG_classes/AnnotSV_ranking', 'Variant_type']
    col.append(f'test0{num}')
    data.columns = col
    return data

def snv_data_arrange(snv_data, num):
    data = pd.read_csv(snv_data)
    col = ['Chr', 'Start', 'End', 'Func.refGene/Location1', 'ExonicFunc.refGene/location2', 'Gene_name', 'AAChange', 'GD_AF_EAS',\
           'drop', 'Variant_ID']
    col.append(f'test0{num}')
    data.columns = col
    data.drop(['drop'], axis=1, inplace=True)
    return data

sv_df = sv_data_arrange(sv_data_list[0], 1)
for i in range(1,len(sv_data_list)):
    df1 = sv_data_arrange(sv_data_list[i], i+1)
    sv_df = pd.merge(sv_df, df1, on=['Chr', 'Start', 'End', 'SV_length', 'SV_type', 'Gene_name', 'Func.refGene/Location1', 'ExonicFunc.refGene/location2',\
                     'Variant_ID', 'GD_AF_EAS', 'GD_POPMAX_AF', 'ACMG_classes/AnnotSV_ranking', 'Variant_type'], how='outer')

snv_df = snv_data_arrange(snv_data_list[0],1)
for i in range(1,len(snv_data_list)):
    df1 = snv_data_arrange((snv_data_list[i]), i+1)
    snv_df = pd.merge(snv_df, df1, on=['Chr', 'Start', 'End', 'Func.refGene/Location1', 'ExonicFunc.refGene/location2', 'Gene_name', 'AAChange', 'GD_AF_EAS',\
                                        'Variant_ID'], how='outer')

data = sv_df.append(snv_df, ignore_index=True)

data_col = ['Chr', 'Start', 'End', 'SV_length', 'SV_type', 'Gene_name', 'Func.refGene/Location1', 'ExonicFunc.refGene/location2', 'Variant_ID',\
            'GD_AF_EAS', 'GD_POPMAX_AF', 'ACMG_classes/AnnotSV_ranking', 'Variant_type', 'AAChange']

for i in range(len(sv_data_list)):
    data_col.append(f'test0{i+1}')
data = data[data_col].fillna('NAN')
sample_df = data[data.columns[14:]]
sample_row = sample_df.shape[0]
Counts = []
for i in range(sample_row):
    t = sample_df.iloc[i].to_list()
    t = [0 if e == '.' else 1 for e in t ]
    count = sum(t)
    Counts.append(count)
data['Counts'] = Counts


# candidate gene filter
candidate = gene[0].to_list()
gene_list = data['Gene_name'].to_list()
candidate_gene_filter = ['+' if gene in candidate else '-' for gene in gene_list]
data['Candidate_gene_filter'] = candidate_gene_filter


# sort columns
data_col = ['Chr', 'Start', 'End', 'SV_length', 'SV_type', 'Gene_name', 'Func.refGene/Location1', 'ExonicFunc.refGene/location2', 'Variant_ID',\
            'GD_AF_EAS', 'GD_POPMAX_AF', 'ACMG_classes/AnnotSV_ranking', 'Variant_type', 'AAChange', 'Candidate_gene_filter']
for i in range(len(sv_data_list)):
    data_col.append(f'test0{i+1}')
data_col.append('Counts')

data = data[data_col]
data.sort_values(by=['Chr', 'Start'], inplace=True)

data.to_csv(args.o, sep='\t', index=False)
