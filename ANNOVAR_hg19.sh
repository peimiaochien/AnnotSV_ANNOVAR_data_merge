#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/staging number
#SBATCH -J hg19_annovar        # Job name
#SBATCH -p ngs48G           # Partition Name 
#SBATCH -c 14               # core preserved 
#SBATCH --mem=46G           # memory used
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL

#hg19 version ANNOVAR with TWB1496 database info, same ANNOVAR in A1_pipeline
#the output result is used for SNV/SV/CNV/repeats mergion due to the python script select the columns by index not by column name. so the order of the columns is critical



tool_dir="/opt/ohpc/Taiwania3/pkg/biology"
logfile=./20220323_run_hg19.log
snv_vcf_input="dragen_v3.9_20220111.hard_filtered.vcf.gz"
output_para=""
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
#set -x


perl ${tool_dir}/ANNOVAR/annovar_20210819/table_annovar.pl ${snv_vcf_input} /staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/ -buildver hg19 -out ${output_para} -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,clinvar_20210123,dbnsfp41a,dbscsnv11,gwava,TWB1496_AF -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f -arg '-splicing 20',,,,,,,,,,,, -nastring . -vcfinput -polish
