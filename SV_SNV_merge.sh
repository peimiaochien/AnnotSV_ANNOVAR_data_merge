#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J merge         # Job name
#SBATCH -p ngs24G           # Partition Name 等同PBS裡面的 -q Queue name
#SBATCH -c 7               
#SBATCH --mem=23G           # 使用的記憶體量 請參考Queue資源設定
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL 




TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x

module load pkg/Anaconda3

data_folder_name="./data"
gene_list="./data/gene.txt"
output="test.txt"

python SV_SNV_merge.py  -dfolder ${data_folder_name}  -g ${gene_list} -o ${output}


