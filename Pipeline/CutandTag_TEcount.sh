#!/bin/sh
#SBATCH --job-name Cut&Tag                
#SBATCH -e TE_count.err
#SBATCH -o TE_count.out
#SBATCH -p your partition
#SBATCH --nodes 1                       
#SBATCH --ntasks-per-node 16          
#SBATCH --cpus-per-task 1                
#SBATCH --time 0-8:00:00               
#SBATCH --mem=20gb


module load htseq
module load samtools


##########################################
################## INPUT #################
##########################################

REFERENCEGFF=dsec_scaffold2_2019.Sech.042623.gff

DIR=/path/to/your/dir/


Species="Sech"

TARGET=("${Species}_CID1" "${Species}_CID2" "${Species}_CID3" "${Species}_H3K9me3")
INPUT="${Species}_IgG"


### Direction
cd ${DIR}

mkdir count_TE


### You need to modify the gff file to be compatible with htseq
awk '{for(i=0;++i<=NF-3;) printf $i"\t" ; print $(NF-2)}'  "$REFERENCEGFF" > ${Species}.htseq.gff


##########################################
################## IgG count #################
##########################################


htseq-count -m intersection-nonempty -f bam -t similarity -i Target Mapping/${INPUT}_bwa_q0.bam ${Species}.htseq.gff > count_TE/${INPUT}.out

##################### Remove the simple repeat
grep -v ')n' count_TE/${INPUT}.out | grep -v '\-rich' | grep -v 'DODECA_SAT' | grep -v '15mer_SAT' > count_TE/${INPUT}_filter.out

##################### Simplify TE names
sed 's/Motif://g' count_TE/${INPUT}_filter.out > count_TE/${INPUT}_filter2.out
rm count_TE/${INPUT}_filter.out

##################### Normalize the count in RPM
##### Number of read mapping for normalization
COUNTq0=$(samtools view -c -F 260 Mapping/${INPUT}_bwa_q0.bam)
NORMq0=$(echo "scale=5;(1000000/$COUNTq0)"| bc -l)
##### Normalization
awk -v a=$NORMq0 '{print $1 "\t" $2*a}' count_TE/${INPUT}_filter2.out | head -n -5 > count_TE/${INPUT}_RPM.out
rm count_TE/${INPUT}_filter2.out


##########################################
############### Target count #################
##########################################


for SAMPLE in "${TARGET[@]}"
do 


##################### Counting the read per TE with htseq
htseq-count -m intersection-nonempty -f bam -t similarity -i Target Mapping/${SAMPLE}_bwa_q0.bam ${Species}.htseq.gff > count_TE/${SAMPLE}.out

##################### Remove the simple repeat
grep -v ')n' count_TE/${SAMPLE}.out | grep -v '\-rich' | grep -v 'DODECA_SAT' | grep -v '15mer_SAT' > count_TE/${SAMPLE}_filter.out
##################### Simplify TE names
sed 's/Motif://g' count_TE/${SAMPLE}_filter.out > count_TE/${SAMPLE}_filter2.out
rm count_TE/${SAMPLE}_filter.out

##################### Normalize the count in RPM

##### Number of read mapping for normalization
COUNTq0=$(samtools view -c -F 260 Mapping/${SAMPLE}_bwa_q0.bam)
NORMq0=$(echo "scale=5;(1000000/$COUNTq0)"| bc -l)

##### Normalization
awk -v a=$NORMq0 '{print $1 "\t" $2*a}' count_TE/${SAMPLE}_filter2.out | head -n -5 > count_TE/${SAMPLE}_RPM.out
rm count_TE/${SAMPLE}_filter2.out

##################### Enrichment : divide by the IgG count
paste count_TE/${Species}_IgG_RPM.out count_TE/${SAMPLE}_RPM.out | awk '{print $1 "\t" $4/($2+1)}' > count_TE/${SAMPLE}_Enrich.out

##################### Got the top 25% more enrich

##### Number correspondig to the top 25% without taking count of the 0 counts
M=$(( $(awk 'NR > 1 {if ($2!=0) print $0}' "count_TE/${SAMPLE}_Enrich.out" | wc -l) / 4 ))

##### Exctract the top 25% from Enrichment file
sort -k2nr count_TE/${SAMPLE}_Enrich.out | head -$M > count_TE/${SAMPLE}_Enrich_Top.out 
sort count_TE/${SAMPLE}_Enrich_Top.out > count_TE/${SAMPLE}_Enrich_Top_sort.out
rm count_TE/${SAMPLE}_Enrich_Top.out 

##### Number correspondig to the top 25% without taking count of the 0 counts
M=$(( $(awk 'NR > 1 {if ($2!=0) print $0}' "count_TE/${SAMPLE}_RPM.out" | wc -l) / 4 ))

##### Exctract the top 25% from RPM file
sort -k2nr count_TE/${SAMPLE}_RPM.out | head -$M > count_TE/${SAMPLE}_RPM_Top.out 
sort count_TE/${SAMPLE}_RPM_Top.out > count_TE/${SAMPLE}_RPM_Top_sort.out
rm  count_TE/${SAMPLE}_RPM_Top.out 


done


##########################################
######## Merge replicate #################
##########################################

### Join all the 3 replicates keeping only the common element

join -1 1 -2 1 count_TE/${Species}_CID1_RPM_Top_sort.out count_TE/${Species}_CID2_RPM_Top_sort.out > count_TE/${Species}_CID_RPM_Top.out_1
join -1 1 -2 1 count_TE/${Species}_CID_RPM_Top.out_1 count_TE/${Species}_CID3_RPM_Top_sort.out > count_TE/${Species}_CID_RPM_Top.out
rm count_TE/${Species}_CID_RPM_Top.out_1

join -1 1 -2 1 count_TE/${Species}_CID1_Enrich_Top_sort.out count_TE/${Species}_CID2_Enrich_Top_sort.out > count_TE/${Species}_CID_Enrich_Top.out_1
join -1 1 -2 1 count_TE/${Species}_CID_Enrich_Top.out_1 count_TE/${Species}_CID3_Enrich_Top_sort.out > count_TE/${Species}_CID_Enrich_Top.out
rm count_TE/${Species}_CID_Enrich_Top.out_1







