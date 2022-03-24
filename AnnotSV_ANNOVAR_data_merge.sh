#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/staging number
#SBATCH -J AnnotSV_ANNOVAR_data_merge         # Job name
#SBATCH -p ngs48G           # Partition Name 
#SBATCH -c 14               # core preserved 
#SBATCH --mem=46G           # memory used
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL



TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${ID}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x


module add pkg/Anaconda3
gene_list="gene.list" 
sv_input="./data/dragen_v3.9_20220323_1.sv.sorted.vcf.annotSV.output.txt"
snv_input="./data/dragen_v3.9_20220323_1.snv.txt"
cnv_input="./data/dragen_v3.9_20220323_1.cnv.sorted.vcf.annotSV.output.txt"
repeats_input="./data/dragen_v3.9_20220323_1.repeats.sorted.vcf.annotSV.output.txt"
output="merging_output.txt"

python  AnnotSV_ANNOVAR_data_merge.py -g ${gene_list} -sv ${sv_input} -snv ${snv_input} -cnv ${cnv_input} -repeats ${repeats_input} -o ${output}
