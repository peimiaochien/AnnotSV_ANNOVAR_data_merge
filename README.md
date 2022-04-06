# AnnotSV_ANNOVAR_data_merge

## Input
1. sample para list
   * one sample should have same para 
3. gene list
   * gene name one is interested in
4. the input files name should include sv/snv/cnv/repeats for the script to distinguish
   * sv/cnv/repeats comes from AnnotSV result 
   * snv comes from ANNOVAR result, use hg19 version ANNOVAR due to TWB1496 database information is needed, the latest ANNOVAR has update in Github already

renew the information in AnnotSV_ANNOVAR_data_merge.sh. 
  sbatch AnnotSV_ANNOVAR_data_merge.sh

## Output



