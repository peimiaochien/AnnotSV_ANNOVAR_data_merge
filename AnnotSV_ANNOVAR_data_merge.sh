#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/staging number
#SBATCH -J AnnotSV_ANNOVAR_data_merge         # Job name
#SBATCH -p ngs24G           # Partition Name 
#SBATCH -c 7               # core preserved 
#SBATCH --mem=23G           # memory used
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
gene_list="./input/gene_list.txt" 
sample_para_list="./input/sample_para_list.txt"
data_folder='./input'
output="./output/merging_output.txt"



# Dont change command bellow
python  AnnotSV_ANNOVAR_data_merge.py -g ${gene_list} -para_list ${sample_para_list} -dfolder ${data_folder} -o ${output}
