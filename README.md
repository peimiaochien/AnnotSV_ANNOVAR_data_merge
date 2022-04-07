# AnnotSV_ANNOVAR_data_merge
The main goal of this script is to merge ANNOVAR and AnnotSV annotation results together, it can be single sample or multiple samples. <br>
Also based on providing gene list to get filter information.

## Input
1. sample para list
   * same sample should have same para 
3. gene list
   * gene name one is interested in
4. the input files name should include sv/snv/cnv/repeats for the script to distinguish
   * sv/cnv/repeats comes from AnnotSV result 
   * snv comes from ANNOVAR result, use hg19 version ANNOVAR due to TWB1496 database information is needed, the latest ANNOVAR has update in Github already

Renew the information in AnnotSV_ANNOVAR_data_merge.sh <br>
`sbatch AnnotSV_ANNOVAR_data_merge.sh` <br>
Base on **Chr, Start, End, Ref, Alt** , the information from AnnotSV and ANNOVAR will merge together.
***
## Output

At the end of the table, there will be columns:
1. **sample column**(column name is given by sample para list) , this column include the genotype of the sample
2. **Counts** how many samples show this variant in this merging data
3. **Candidate_gene_filter** if this variant take place on the gene in candidate gene list than it will mark as '+', else mark as '-'


