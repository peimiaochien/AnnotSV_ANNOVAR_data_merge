#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/staging number
#SBATCH -J AnnotSV_ANNOVAR_data_merge         # Job name
#SBATCH -p ngs96G           # Partition Name 
#SBATCH -c 28               # core preserved 
#SBATCH --mem=92G           # memory used
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL


# This script goal is to merge ANNOVAR and AnnotSV results together from one or multiple samples
# Also use gene list to make filter



TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${ID}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x

# update your files path here
module add pkg/Anaconda3
gene_list="./gene_list.txt" 
sample_para_list="./sample_para_list.txt"
data_folder='./input'
output="./merging_output.txt"
# reference_version only hg38 or hg19
reference_version='hg38'


# Dont change command bellow
python  AnnotSV_ANNOVAR_data_merge.py -g ${gene_list} -para_list ${sample_para_list} -dfolder ${data_folder} -o ${output} -ref ${reference_version}
