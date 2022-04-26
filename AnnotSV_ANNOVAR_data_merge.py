import pandas as pd
import sys
import glob
import argparse
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
parser = argparse.ArgumentParser()
parser.add_argument("-dfolder", help="put all data you want to merge in the folder, same sample merging files must have same para and in  use SV/SNV/repeat/CNV to annotate file type in file name")
parser.add_argument("-g", help="gene list in csv formate")
parser.add_argument("-o", help="output file name")
parser.add_argument("-para_list", help='sample_para_list')
parser.add_argument('-ref', help='reference genome version ANNOVAR use, hg19 or hg38')
args = parser.parse_args()


candidate_gene_list = pd.read_csv(args.g, sep='\t', header=None)[0].to_list()

def annovar_data_arrange(annovar_data, sample_para, ref_version):
    if ref_version == 'hg19':
        data = pd.read_csv(annovar_data, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 8, 9, 17, 21, 22, 27, 29, 38, 39, 40, 41, 46, 48, 118, 119, 129, 130, 131, 132])
        data.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
           'ExonicFunc.refGene', 'AAChange.refGene', 'Gene_name', 'genome_AF',
           'genome_AF_popmax', 'genome_AF_sas', 'genome_AF_eas', 'avsnp150',
           'TaiwanBiobank-official_Illumina1000-AF', 'exome_AF', 'exome_popmax',
           'exome_AF_sas', 'exome_AF_eas', 'TWB1496_AF', 'TWB1496_QC', 'Otherinfo10',
           'Otherinfo11', 'Otherinfo12', 'Otherinfo13']
    elif ref_version == 'hg38':
        data = pd.read_csv(annovar_data, sep='\t', usecols=[0, 1, 2, 3, 4, 5, 6, 8, 9, 17, 21, 29, 33, 34, 35, 36, 41, 43, 119, 120, 121, 122])
        data.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene_name', 'ExonicFunc.refGene', 'AAChange.refGene', 'Gene.ensGene', 'genome_AF', 'genome_AF_eas', 'genome_AF_sas',
                        'avsnp150', 'exome_AF', 'exome_AF_popmax', 'exome_AF_sas', 'exome_AF_eas', 'Otherinfo10', 'Otherinfo11', 'Otherinfo12', 'Otherinfo13']
    else:
        print('Only provide hg19 and hg38 ANNOVAR annotation results in merging')
    data['Chr'] = data['Chr'].astype('object')
    data[sample_para] = [ variant.split(',')[0].split(':')[0] for variant in data['Otherinfo13'].to_list()]
    return data


def annotsv_data_arrange(annotsv_data, sample_para):
    data = pd.read_csv(annotsv_data, sep='\t', usecols=[1,2,3,4,5,7,8,9,10, 11, 12, 13, 14, 15, 16, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95], low_memory=False)
    data.columns = ['Chr', 'Start', 'End', 'SV_length', 'SV_type', 'ID', 'Ref', 'Alt', 'QUAL',
       'FILTER', 'INFO', 'FORMAT', 'sample', 'Annotation_mode', 'Gene_name',
       'Tx', 'Tx_start', 'Tx_end', 'Overlapped_tx_length',
       'Overlapped_CDS_length', 'Overlapped_CDS_percent', 'Frameshift',
       'Exon_count', 'Location', 'Location2', 'Dist_nearest_SS',
       'Nearest_SS_type', 'Intersect_start', 'Intersect_end', 'RE_gene',
       'P_gain_phen', 'P_gain_hpo', 'P_gain_source', 'P_gain_coord',
       'P_loss_phen', 'P_loss_hpo', 'P_loss_source', 'P_loss_coord',
       'P_ins_phen', 'P_ins_hpo', 'P_ins_source', 'P_ins_coord',
       'P_snvindel_nb', 'P_snvindel_phen', 'B_gain_source', 'B_gain_coord',
       'B_loss_source', 'B_loss_coord', 'B_ins_source', 'B_ins_coord',
       'B_inv_source', 'B_inv_coord', 'TAD_coordinate', 'ENCODE_experiment',
       'GC_content_left', 'GC_content_right', 'Repeat_coord_left',
       'Repeat_type_left', 'Repeat_coord_right', 'Repeat_type_right',
       'Gap_left', 'Gap_right', 'SegDup_left', 'SegDup_right',
       'ENCODE_blacklist_left', 'ENCODE_blacklist_characteristics_left',
       'ENCODE_blacklist_right', 'ENCODE_blacklist_characteristics_right',
       'ACMG', 'HI', 'TS', 'DDD_HI_percent', 'DDD_status', 'DDD_mode',
       'DDD_consequence', 'DDD_disease', 'DDD_pmid', 'ExAC_delZ', 'ExAC_dupZ',
       'ExAC_cnvZ', 'ExAC_synZ', 'ExAC_misZ', 'OMIM_ID', 'OMIM_phenotype',
       'OMIM_inheritance', 'OMIM_morbid', 'OMIM_morbid_candidate', 'LOEUF_bin',
       'GnomAD_pLI', 'ExAC_pLI', 'AnnotSV_ranking_score',
       'AnnotSV_ranking_criteria', 'ACMG_class']

    data = data[data['Annotation_mode'] == 'split']
    data['Chr'] = [f'chr{num}' for num in data['Chr'].to_list()]
    data[sample_para] = [sample.split(':')[0] for sample in data['sample'].to_list()]
    data.drop(['sample'], axis=1, inplace=True)
    return data



def annotsv(sample_para):
    sv_data = glob.glob(f'{args.dfolder}/{sample_para}*sv*')[0]
    cnv_data = glob.glob(f'{args.dfolder}/{sample_para}*cnv*')[0]
    repeats_data = glob.glob(f'{args.dfolder}/{sample_para}*cnv*')[0]
    sv_df = annotsv_data_arrange(sv_data, sample_para)
    cnv_df = annotsv_data_arrange(cnv_data, sample_para)
    repeats_df = annotsv_data_arrange(repeats_data, sample_para)
    return sv_df, cnv_df, repeats_df

def annovar(sample_para, ref_version):
    snv_data = glob.glob(f'{args.dfolder}/{sample_para}*snv*')[0]
    snv_df = annovar_data_arrange(snv_data, sample_para, ref_version)
    return snv_df

#sample import
total_df = pd.DataFrame()
sample_list = pd.read_csv(args.para_list, sep='\t', header=None)[0].to_list()
for sample_para in sample_list:
    sv_df, cnv_df, repeats_df = annotsv(sample_para)
    snv_df = annovar(sample_para, args.ref)
    total_df = total_df.append([snv_df, sv_df, cnv_df, repeats_df], ignore_index=True)
total_df.fillna(value='.', inplace=True)

col_list = total_df.columns.to_list()[5:]
for sample_para in sample_list:
    col_list.remove(sample_para)

total_df = total_df.groupby(['Chr', 'Start', 'End', 'Ref', 'Alt'])

final_data = total_df.size().to_frame(name='count').reset_index().drop(['count'], axis=1)
# other than sample
for col in col_list:
    final_data[col] = [data[0] for data in total_df[col].agg([(col, lambda x: list(set(x)))])[col].tolist()]

# data
for sample_para in sample_list:
    final_data[sample_para] = [ data[1] if len(data) > 1 else data[0] for data in total_df[sample_para].agg([(sample_para, lambda x:list(set(x)))])[sample_para].tolist()]

# count
Counts = []
for i in range(len(final_data)):
    count=0
    for sample in sample_list:
        if (final_data[sample].iloc[i] != 'NaN') and (final_data[sample].iloc[i] != '.') and (final_data[sample].iloc[i] != '0/0'):
            count += 1
    Counts.append(count)
final_data['Counts'] = Counts

# candidate gene filter

gene_list = final_data['Gene_name'].to_list()
candidate_gene_filter = ['+' if gene in candidate_gene_list else '-' for gene in gene_list]
final_data['Candidate_gene_filter'] = candidate_gene_filter

final_data.sort_values(by=['Chr', 'Start'], inplace=True)
final_data.to_csv(args.o, sep='\t', index=False)


