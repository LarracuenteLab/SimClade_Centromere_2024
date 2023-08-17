#!/bin/sh
#SBATCH --job-name blast_sat             
#SBATCH -e blast_sat.err
#SBATCH -o blast_sat.out
#SBATCH -p debug
#SBATCH -c 2
#SBATCH -n 1                      
#SBATCH --time 1:00:00                
#SBATCH --mem=10gb


module load blast


### Direction
DIR="/path/to/your/dir/"
cd ${DIR}


##########################################
################## INPUT #################
##########################################

input="365bp.fasta"

REF=dmel_scaffold2_plus0310.fasta
Species="Dmel"

##########################################
################## blast #################
##########################################

makeblastdb -in ${REF} -input_type fasta -dbtype nucl -out ${Species}

blastn -task blastn -db ${Species} -query ${input}  -outfmt 6 -out ${Species}_${input}.out 
