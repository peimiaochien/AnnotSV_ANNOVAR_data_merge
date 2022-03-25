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
gene_list="gene_list.txt" 
sample_para_list="sample_para_list.txt"
data_folder='./data'
output="merging_output.txt"




python  AnnotSV_ANNOVAR_data_merge.py -g ${gene_list} -para_list ${sample_para_list} -dfolder ${data_folder} -o ${output}
