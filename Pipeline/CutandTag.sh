#!/bin/sh
#SBATCH --job-name Cut&Tag                
#SBATCH -p your partition
#SBATCH -e Cut&Tag.err
#SBATCH -o Cut&Tag.out
#SBATCH --nodes 1                       
#SBATCH --ntasks-per-node 20          
#SBATCH --cpus-per-task 1                
#SBATCH --time 0-48:00:00                
#SBATCH --mem=120gb


##########################################
################## INPUT #################
##########################################
REF=dmel_scaffold2_plus0310.fasta

DIR=/path/to/your/dir/


Species="N25"

TARGET=("${Species}_CID1" "${Species}_CID2" "${Species}_CID3" "${Species}_H3K9me3" "${Species}_IgG")


cd ${DIR}


#####################
######Index #########
######################
module load bwa
bwa index ${REF} 


##########################################
#################  PIPELINE ##############
##########################################

for SAMPLE in "${TARGET[@]}"
do


module purge
module load bwa
module load samtools
module load bedtools
module load trimgalore/0.4.4 cutadapt/b1 htseq/0.9.1


#################################
########## Trimming ##########
#################################
mkdir fastq_Trim

if [ -f "fastq_Trim/${SAMPLE}_R1_val_1.fq.gz" ]; then
	    echo "reads already trimmed"
else	
		echo "trimming reads"

trim_galore --paired --nextera --length 75 --phred33 --no_report_file --fastqc -o fastq_Trim fastq/${SAMPLE}_R1.fastq.gz fastq/${SAMPLE}_R2.fastq.gz

fi

N=$(( $(zcat "fastq_Trim/${SAMPLE}_R1_val_1.fq.gz" | wc -l) /4 ))
M=$(( $(zcat "fastq_Trim/${SAMPLE}_R2_val_2.fq.gz" | wc -l) /4 ))
	
	echo "${SAMPLE} number of raw reads: R1: $N; R2: $M"

#################################
########## Mapping ##############
#################################

if [ -f "Mapping/${SAMPLE}_bwa_sorted.bam" ]; then
	echo "reads already mapped"
else
	echo "Mapping"
	

########### Mapping ########### 
 
mkdir Mapping

bwa mem -t 20 ${REF} \
fastq_Trim/${SAMPLE}_R1_val_1.fq.gz  \
fastq_Trim/${SAMPLE}_R2_val_2.fq.gz \
> Mapping/${SAMPLE}_bwa.sam


samtools view -bSh Mapping/${SAMPLE}_bwa.sam  >  Mapping/${SAMPLE}_bwa.bam
rm Mapping/${SAMPLE}_bwa.sam 

samtools sort -@ 20 Mapping/${SAMPLE}_bwa.bam > Mapping/${SAMPLE}_bwa_sorted.bam
rm Mapping/${SAMPLE}_bwa.bam

samtools index Mapping/${SAMPLE}_bwa_sorted.bam

samtools flagstat  Mapping/${SAMPLE}_bwa_sorted.bam > Mapping/${SAMPLE}_bwa_sorted.bam.out

### Mapping quality
samtools view Mapping/${SAMPLE}_bwa_sorted.bam | awk -F "\t" '{print $5}' > Mapping/${SAMPLE}_bwa_sorted.txt

fi

N=$(( $(samtools view "Mapping/${SAMPLE}_bwa_sorted.bam" | wc -l) /2 ))
echo "${SAMPLE} got mapped using bwa total pairs mapped: $N"


###########################################################################
############################## Mark Duplicate ####################################
###########################################################################

if [ -f "Mapping/${SAMPLE}_bwa_noDup.bam" ]; then
	echo "Duplicate alread removed"
else
	echo "Remove duplicate"
	
java -Xmx10g -jar /software/picard/2.12.0/picard.jar MarkDuplicates \
          REMOVE_DUPLICATES=TRUE \
          I=Mapping/${SAMPLE}_bwa_sorted.bam \
          O=Mapping/${SAMPLE}_bwa_noDup.bam \
          M=Mapping/${SAMPLE}_bwa_noDup.txt

fi 


###########################################################################
######################### Quality filtering ################################
###########################################################################


if [ -f "Mapping/${SAMPLE}_bwa_q30.bam" ]; then
	echo "reads filter"
else
	echo "Quality Filtering"

### Filtering : keeping multi-mapping reads
samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048 Mapping/${SAMPLE}_bwa_noDup.bam > Mapping/${SAMPLE}_bwa_q0.bam
samtools index Mapping/${SAMPLE}_bwa_q0.bam

N=$(( $(samtools view "Mapping/${SAMPLE}_bwa_q0.bam" | wc -l) /2 ))
echo "${SAMPLE}  total pairs multi apped after filtering : $N"



### Filtering : keeping uniquely mapping read
samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 2048 -q 30 Mapping/${SAMPLE}_bwa_noDup.bam > Mapping/${SAMPLE}_bwa_q30.bam
samtools index  Mapping/${SAMPLE}_bwa_q30.bam

N=$(( $(samtools view "Mapping/${SAMPLE}_bwa_q30.bam" | wc -l) /2 ))
echo "${SAMPLE}  total pairs uniquely mapped after filtering : $N"


rm Mapping/${SAMPLE}_bwa_noDup.bam
fi

###########################################################################
######################### Coverage ################################
###########################################################################



if [ -f "Coverage/${SAMPLE}_bwa_q0_RPM.bw" ]; then
	echo "Coverage : done "
else
	echo "Coverage"

mkdir Coverage

##### Number of read mapping for normalization
COUNTq0=$(samtools view -c -F 260 Mapping/${SAMPLE}_bwa_q0.bam)
NORMq0=$(echo "scale=5;(1000000/$COUNTq0)"| bc -l)

###### Coverage
bamCoverage -b Mapping/${SAMPLE}_bwa_q0.bam -o Coverage/${SAMPLE}_bwa_q0_RPM.bw --scaleFactor ${NORMq0} -bs 1 --extendReads

##### Number of read mapping for normalization
COUNTq30=$(samtools view -c -F 260 Mapping/${SAMPLE}_bwa_q30.bam)
NORMq30=$(echo "scale=5;(1000000/$COUNTq30)"| bc -l)

###### Coverage
bamCoverage -b Mapping/${SAMPLE}_bwa_q30.bam -o Coverage/${SAMPLE}_bwa_q30_RPM.bw --scaleFactor ${NORMq30} -bs 1 --extendReads

fi

done


