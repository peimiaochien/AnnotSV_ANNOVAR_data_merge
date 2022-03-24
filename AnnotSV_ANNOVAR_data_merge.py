import pandas as pd
import sys
import glob
import argparse
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
parser = argparse.ArgumentParser()
parser.add_argument("-sv")
parser.add_argument("-snv")
parser.add_argument("-cnv")
parser.add_argument("-repeats")
parser.add_argument("-g", help="gene list in csv formate")
parser.add_argument("-o", help="output file name")
args = parser.parse_args()


# data import
sv_data_list = args.sv
snv_data_list = args.snv
repeat_data_list = args.repeats
cnv_data_list = args.cnv

gene = pd.read_csv(args.g, sep='\t', header=None)
sample_list = []

# annovar multianno.txt
def annovar_data_arrange(annovar_data):
    data = pd.read_csv(annovar_data, sep='\t',
                       usecols=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
                                'ExonicFunc.refGene', 'AAChange.refGene', 'Gene.ensGene', 'AF',
                                'AF_popmax', 'AF_sas', 'AF_eas', 'avsnp150',
                                'TaiwanBiobank-official_Illumina1000-AF', 'AF.1', 'AF_popmax.1',
                                'AF_sas.1', 'AF_eas.1', 'TWB1496_AF', 'TWB1496_QC', 'Otherinfo10',
                                'Otherinfo11', 'Otherinfo12', 'Otherinfo13'])
    data.columns = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
       'ExonicFunc.refGene', 'AAChange.refGene', 'Gene.ensGene', 'genome_AF',
       'genome_AF_popmax', 'genome_AF_sas', 'genome_AF_eas', 'avsnp150',
       'TaiwanBiobank-official_Illumina1000-AF', 'exome_AF', 'exome_popmax',
       'exome_AF_sas', 'exome_AF_eas', 'TWB1496_AF', 'TWB1496_QC', 'Otherinfo10',
       'Otherinfo11', 'Otherinfo12', 'Otherinfo13']
    data['Chr'] = data['Chr'].astype('object')
    snv_para = annovar_data.split('/')[-1]
    sample_list.append(snv_para)
    data[snv_para] = [ variant.split(',')[0].split(':')[0] for variant in data['Otherinfo13'].to_list()]
    data.drop(['Otherinfo13'], axis=1, inplace=True)
    return data

# re-asign columns 20220321
def annotsv_data_arrange(annotsv_data):
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
    data['Chr'] = [f"chr{num}" for num in data['Chr'].tolist()]
    sv_para = annotsv_data.split('/')[-1]
    sample_list.append(sv_para)
    data[sv_para] = [sample.split(':')[0] for sample in data['sample'].to_list()]
    data.drop(['sample'], axis=1, inplace=True)
    return data

sv_df = annotsv_data_arrange(sv_data_list)
repeat_df = annotsv_data_arrange(repeat_data_list)
cnv_df = annotsv_data_arrange(cnv_data_list)
snv_df = annovar_data_arrange(snv_data_list)



# append --> concat method
data = sv_df.merge(snv_df, on=['Chr', 'Start', 'End','Ref', 'Alt'], how='outer')
data = data.merge(cnv_df, on=['Chr', 'Start', 'End','Ref', 'Alt'], how='outer')
data = data.merge(repeat_df, on=['Chr', 'Start', 'End','Ref', 'Alt'], how='outer')
data.fillna(value='NaN', inplace=True)

### 20220322
Counts = []
for i in range(len(data)):
    count=0
    for sample in sample_list:
        if (data[sample].iloc[i] != 'NaN') and (data[sample].iloc[i] != '.'):
            count += 1
    Counts.append(count)

data['Counts'] = Counts


# candidate gene filter
candidate = gene[0].to_list()
gene_list = data['Gene_name'].to_list()
candidate_gene_filter = ['+' if gene in candidate else '-' for gene in gene_list]
data['Candidate_gene_filter'] = candidate_gene_filter


data.sort_values(by=['Chr', 'Start'], inplace=True)
data.to_csv(args.o, sep='\t', index=False)
