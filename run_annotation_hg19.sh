#!/usr/bin/bash
#SBATCH -A MST109178        # Account name/project number
#SBATCH -J annovar         # Job name
#SBATCH -p ngs96G           # Partition Name 
#SBATCH -c 28               # core preserved 
#SBATCH --mem=92G           # memory used
#SBATCH -o out.log          # Path to the standard output file 
#SBATCH -e err.log          # Path to the standard error ouput file
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL 
### Please define the following variables
INPUT="/staging/reserve/paylong_ntu/AI_SHARE/GitHub/ANNOTATION/TestData/DF_Fid.vqsr_SNP_INDEL.hc.recaled.vcf.gz"
wkdir="/staging/reserve/paylong_ntu/AI_SHARE/GitHub/ANNOTATION/ANNOVAR/"
para="github_ANNOVAR_hg19"

### DO NOT CHANGE
REF="/staging/reserve/paylong_ntu/AI_SHARE/reference/GATK_bundle/2.8/hg19/ucsc.hg19.fasta"
ANNOVAR="/opt/ohpc/Taiwania3/pkg/biology/ANNOVAR/annovar_20210819/table_annovar.pl"
humandb="/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb"
QC_ADDING_SH="TWB1496_QC_adding.sh"
hg19_TWB1496_QC='/staging/reserve/paylong_ntu/AI_SHARE/reference/annovar_2016Feb01/humandb/hg19_TWB1496_QC.txt'
cd ${wkdir}
mkdir -p ${wkdir}
cd ${wkdir}
TIME=`date +%Y%m%d%H%M`
logfile=./${TIME}_${para}_run.log
exec 3<&1 4<&2
exec >$logfile 2>&1
set -euo pipefail
set -x

/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools norm -m- $INPUT -O z -o ${wkdir}/${para}.decom_hg19.vcf.gz
/opt/ohpc/Taiwania3/pkg/biology/BCFtools/bcftools_v1.13/bin/bcftools norm -f $REF ${wkdir}/${para}.decom_hg19.vcf.gz -O z -o ${wkdir}/${para}.decom_hg19.norm.vcf.gz

perl $ANNOVAR ${wkdir}/${para}.decom_hg19.norm.vcf.gz ${humandb} -buildver hg19 -out ${para} -remove -protocol refGene,cytoBand,knownGene,ensGene,gnomad211_genome,avsnp150,TaiwanBiobank-official,gnomad211_exome,clinvar_20210123,dbnsfp41a,dbscsnv11,gwava,TWB1496_AF -operation gx,r,gx,gx,f,f,f,f,f,f,f,f,f -arg '-splicing 20',,,,,,,,,,,,,,,,,  -nastring . -vcfinput -polish
head -n 1 ${para}.hg19_multianno.txt > ${para}.filtered_annotation.txt
#grep -P "\texonic\t" ${para}.hg19_multianno.txt | grep -P -v "\tsynonymous" >> ${para}.filtered_annotation.txt
grep -e exonic -e splicing ${para}.hg19_multianno.txt | grep -P -v "\tsynonymous" | grep -P -v "\tncRNA_exonic\t" >> ${para}.filtered_annotation.txt

#Adding hg19_TWB1496 qc
sh ${QC_ADDING_SH} -a ${para}.filtered_annotation.txt -q ${hg19_TWB1496_QC} -o ${para}.filtered_annotation_QC.txt



#done</work2/u1067478/Annotation/CGMH_AML/CGMH_AML_NameList.txt
rm ${para}.avinput ${para}.decom_hg19.norm.vcf.gz ${para}.decom_hg19.vcf.gz ${para}.hg19_multianno.vcf
