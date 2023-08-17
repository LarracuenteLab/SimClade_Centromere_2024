#!/bin/sh
#SBATCH --job-name PeakCalling                
#SBATCH -p your partition
#SBATCH --nodes 1                       
#SBATCH --ntasks-per-node 20          
#SBATCH --cpus-per-task 1                
#SBATCH --time 0-8:00:00                
#SBATCH --mem=20gb


cd /path/to/your/dir/

module load macs/2017-10-26

##########################################
################## IMPUT #################
##########################################


Species="XD1"

TARGET=("${Species}_CID1" "${Species}_CID2" "${Species}_CID3" "${Species}_H3K9me3")

IgG="${Species}_IgG"

###########################################################################
#################### Peak Calling #############################
###########################################################################

mkdir Peaks

for SAMPLE in "${TARGET[@]}"
do


if [ -f "Peaks/${SAMPLE}_bwa_q30_macs2_peaks.narrowPeak" ]; then
	echo "Peaks calling alread done"
else
	echo "Peak Calling"
	
macs2 callpeak -t Mapping/${SAMPLE}_bwa_q30.bam -c Mapping/${IgG}_bwa_q30.bam -f BAMPE -n ${SAMPLE}_bwa_q30_macs2 -g dm -q 0.01 -B --outdir Peaks --call-summits
macs2 callpeak -t Mapping/${SAMPLE}_bwa_q0.bam -c Mapping/${IgG}_bwa_q0.bam -f BAMPE -n ${SAMPLE}_bwa_q0_macs2 -g dm -q 0.01 -B --outdir Peaks --call-summits

fi 

done


###########################################################################
#################### IDR #############################
###########################################################################
module purge
#module load python3/3.7.7
module load idr/2.0.3

cd Peaks

#Sort peak by -log10(p-value)
sort -k8,8nr ${Species}_CID1_bwa_q30_macs2_peaks.narrowPeak > ${Species}_CID1_bwa_q30_macs2_sort.peaks.narrowPeak
sort -k8,8nr ${Species}_CID2_bwa_q30_macs2_peaks.narrowPeak > ${Species}_CID2_bwa_q30_macs2_sort.peaks.narrowPeak
sort -k8,8nr ${Species}_CID3_bwa_q30_macs2_peaks.narrowPeak > ${Species}_CID3_bwa_q30_macs2_sort.peaks.narrowPeak

idr --samples ${Species}_CID1_bwa_q30_macs2_sort.peaks.narrowPeak ${Species}_CID2_bwa_q30_macs2_sort.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ${Species}.1.q30-idr \
--plot \
--log-output-file ${Species}.1.q30.idr.log

idr --samples ${Species}_CID1_bwa_q30_macs2_sort.peaks.narrowPeak ${Species}_CID3_bwa_q30_macs2_sort.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ${Species}.2.q30-idr \
--plot \
--log-output-file ${Species}.2.q30.idr.log

idr --samples ${Species}.1.q30-idr ${Species}.2.q30-idr \
--input-file-type narrowPeak \
--rank p.value \
--output-file ${Species}.q30-idr \
--plot \
--log-output-file ${Species}.q30.idr.log

### taking only the peak with IDR >540
awk '{if($5 >= 540) print $0}' ${Species}.q30-idr > ${Species}.q30.Top-idr


sort -k8,8nr ${Species}_CID1_bwa_q0_macs2_peaks.narrowPeak > ${Species}_CID1_bwa_q0_macs2_sort.peaks.narrowPeak
sort -k8,8nr ${Species}_CID2_bwa_q0_macs2_peaks.narrowPeak > ${Species}_CID2_bwa_q0_macs2_sort.peaks.narrowPeak
sort -k8,8nr ${Species}_CID3_bwa_q0_macs2_peaks.narrowPeak > ${Species}_CID3_bwa_q0_macs2_sort.peaks.narrowPeak

idr --samples ${Species}_CID1_bwa_q0_macs2_sort.peaks.narrowPeak ${Species}_CID2_bwa_q0_macs2_sort.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ${Species}.1.q0-idr \
--plot \
--log-output-file ${Species}.1.q0.idr.log

idr --samples ${Species}_CID1_bwa_q0_macs2_sort.peaks.narrowPeak ${Species}_CID3_bwa_q0_macs2_sort.peaks.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file ${Species}.2.q0-idr \
--plot \
--log-output-file ${Species}.2.q0.idr.log

idr --samples ${Species}.1.q0-idr ${Species}.2.q0-idr \
--input-file-type narrowPeak \
--rank p.value \
--output-file ${Species}.q0-idr \
--plot \
--log-output-file ${Species}.q0.idr.log


### taking only the peak with IDR >540
awk '{if($5 >= 540) print $0}' ${Species}.q0-idr > ${Species}.q0.Top-idr
