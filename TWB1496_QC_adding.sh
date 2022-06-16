#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/staging number
#SBATCH -J A1_panel         # Job name
#SBATCH -p ngs24G          # Partition Name 
#SBATCH -c 7               # core preserved 
#SBATCH --mem=23G           # memory used
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL



# get arguments
while getopts a:q:o: flag
do
    case "${flag}" in
    a) an_file=${OPTARG};;
    q) qc_file=${OPTARG};;
    o) output=${OPTARG};;
    esac
done


qc_row=$(awk '-F\t' '{for (i=1;i<=NF;++i){if ($i ~ "TWB1496_QC")print i}}' ${an_file} )


fun1(){
 l=$1
 awk '-F\t' 'NR==l{print}' l=${l} ${an_file} > "an${l}.tmp"
 if [[ "${an_qc}" == "." ]];
 then
     awk '-F\t' 'st>=$2&&$3>=ed{print$0}' st=${st} ed=${ed} ${tmp_file} > "qc${l}.tmp"
     if [[ $(wc -l qc${l}.tmp| awk '-F ' '{print$1}') == 0 ]];
     then
         awk '-F\t' 'st>=$2&&$3>=st{print$0}' st=${st} ed=${ed} ${tmp_file} > "qc${l}.tmp"
         awk '-F\t' 'ed>=$2&&$3>=ed {print$0}' st=${st} ed=${ed} ${tmp_file} >> "qc${l}.tmp"
     fi

## same chr start end 
    if [[ $(wc -l qc${l}.tmp| awk '-F ' '{print$1}') == 0 ]];
    then
        qc="ND"
    else
## qc is overlap by the data
	qc="+"$(awk '-F\t' '{print$4}' qc${l}.tmp| sed ':a;N;$!ba;s/\n/,/g')
    fi  
    awk '-F\t' -v OFS="\t" '{$qc_row=qc}1' qc=${qc} qc_row=${qc_row}  an${l}.tmp >> "${output}"
 else
     cat an${l}.tmp >> "${output}"
 fi

rm qc${l}.tmp an${l}.tmp
}

i=1



while IFS=$'\t' read -r $(head -n 1 ${an_file}|sed 's/\./_/g;s/-/_/g;s/+/_/g'); 
do
 chr='chr'${Chr}
 st=${Start}
 ed=${End}
 an_qc=${TWB1496_QC}
 # create own qc file#
 export tmp_file=${chr}.tmp
 if [[ -f "${tmp_file}" ]];
 then 
     echo ${tmp_file}
 else
     awk  '-F\t' '$1==chr{print$0}' chr=${chr} ${qc_file} > ${tmp_file}
 fi

 fun1 ${i} &

 (( i++ ))


done< ${an_file} 

wait

echo "all done"
