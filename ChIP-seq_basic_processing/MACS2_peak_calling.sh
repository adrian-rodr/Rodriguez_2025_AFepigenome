#!/bin/bash
#$ -cwd           # Set the working directory for the job to the current directory
#$ -j y           # Join stdout and stderr
#$ -pe smp 6       # Request 1 CPU core
#$ -l h_rt=24:0:0  # Request 1 hour runtime
#$ -l h_vmem=20G   # Request 1GB RAM / core, i.e. 1GB total


cd /path/to/lab/directory/ #Set working directory to path containing preprocessed BAM files from BAM_files_preprocessing.pl

#load python and macs2 (via virtual environment)
#module load python/3.6.3
module load python
#virtualenv --include-lib macs2
virtualenv macs2
source macs2/bin/activate
pip install numpy
pip install macs2


# Remove duplicates with Picard in BAM files from libraries displaying high levels of duplication

module load samtools
samtools view -bF 4 vl0214.bam > vl0214_filtered.bam
##MarkDuplicates:
module load java/1.8.0_121-oracle
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0214_filtered.bam -O vl0214_filtered_dupremoved.bam -M vl0214_marked_dup_metrics.txt

samtools view -bF 4 vl0217.bam > vl0217_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0217_filtered.bam -O vl0217_filtered_dupremoved.bam -M vl0217_marked_dup_metrics.txt

samtools view -bF 4 vl0261.bam > vl0261_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0261_filtered.bam -O vl0261_filtered_dupremoved.bam -M vl0261_marked_dup_metrics.txt

samtools view -bF 4 vl0323.bam > vl0323_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0323_filtered.bam -O vl0323_filtered_dupremoved.bam -M vl0323_marked_dup_metrics.txt

samtools view -bF 4 vl0326.bam > vl0326_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0326_filtered.bam -O vl0326_filtered_dupremoved.bam -M vl0326_marked_dup_metrics.txt

samtools view -bF 4 vl0327.bam > vl0327_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0327_filtered.bam -O vl0327_filtered_dupremoved.bam -M vl0327_marked_dup_metrics.txt

samtools view -bF 4 vl0328.bam > vl0328_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0328_filtered.bam -O vl0328_filtered_dupremoved.bam -M vl0328_marked_dup_metrics.txt

samtools view -bF 4 vl0333.bam > vl0333_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0333_filtered.bam -O vl0333_filtered_dupremoved.bam -M vl0333_marked_dup_metrics.txt

samtools view -bF 4 vl0335.bam > vl0335_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0335_filtered.bam -O vl0335_filtered_dupremoved.bam -M vl0335_marked_dup_metrics.txt

samtools view -bF 4 vl0339.bam > vl0339_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0339_filtered.bam -O vl0339_filtered_dupremoved.bam -M vl0339_marked_dup_metrics.txt

samtools view -bF 4 vl0340.bam > vl0340_filtered.bam
java -jar picard.jar MarkDuplicates -CREATE_INDEX true -REMOVE_DUPLICATES true -I vl0340_filtered.bam -O vl0340_filtered_dupremoved.bam -M vl0340_marked_dup_metrics.txt




# Declare variables for samples in AF cohort
vl0026_12b_SR_LA_H3K27ac="vl0026_20M.bam"
vl0027_12b_SR_RA_H3K27ac="vl0027_20M.bam"
vl0028_25_AF_LA_H3K27ac="vl0028_20M.bam"
vl0029_25_AF_RA_H3K27ac="vl0029_20M.bam"
vl0031_27_SR_RA_H3K27ac="vl0031_20M.bam"
vl0032_28_AF_LA_H3K27ac="vl0032_20M.bam"
vl0033_28_AF_RA_H3K27ac="vl0033_20M.bam"
vl0034_12b_SR_LA_Input="vl0034_20M.bam"
vl0035_12b_SR_RA_Input="vl0035_20M.bam"
vl0036_25_AF_LA_Input="vl0036_20M.bam"
vl0037_25_AF_RA_Input="vl0037_20M.bam"
vl0039_27_SR_RA_Input="vl0039_20M.bam"
vl0040_28_AF_LA_Input="vl0040_20M.bam"
vl0041_28_AF_RA_Input="vl0041_20M.bam"
vl0212_47_AF_LA_H3K27ac="vl0212_20M.bam"
vl0213_45_SR_LA_H3K27ac="vl0213_20M.bam"
vl0214_2231_AF_RA_H3K27ac="vl0214_filtered_dupremoved.bam"
vl0215_31_AF_LA_H3K27ac="vl0215_20M.bam"
vl0217_20_SR_LA_H3K27ac="vl0217_filtered_dupremoved.bam"
vl0218_33_SR_RA_H3K27ac="vl0218_20M.bam"
vl0228_47_AF_LA_Input="vl0228_20M.bam"
vl0229_45_SR_LA_Input="vl0229_20M.bam"
vl0231_31_AF_LA_Input="vl0231_20M.bam"
vl0234_33_SR_RA_Input="vl0234_20M.bam"
vl0256_2231_AF_LA_H3K27ac="vl0256_20M.bam"
vl0257_313_AF_LA_H3K27ac="vl0257_20M.bam"
vl0258_31_AF_RA_H3K27ac="vl0258_20M.bam"
vl0260_313_AF_RA_H3K27ac="vl0260_20M.bam"
vl0261_20_SR_RA_H3K27ac="vl0261_filtered_dupremoved.bam"
vl0262_2231_AF_LA_Input="vl0262_20M.bam"
vl0263_313_AF_LA_Input="vl0263_20M.bam"
vl0264_31_AF_RA_Input="vl0264_20M.bam"
vl0266_313_AF_RA_Input="vl0266_20M.bam"
vl0267_20_SR_RA_Input="vl0267_filtered_dupremoved.bam"
vl0326_18_AF_RA_H3K27ac="vl0326_filtered_dupremoved.bam"
vl0328_27_SR_LA_H3K27ac="vl0328_filtered_dupremoved.bam"
vl0333_18_AF_RA_Input="vl0333_filtered_dupremoved.bam"
vl0335_27_SR_LA_Input="vl0335_filtered_dupremoved.bam"
vl0339_2231_AF_RA_Input="vl0339_filtered_dupremoved.bam"
vl0340_20_SR_LA_Input="vl0340_filtered_dupremoved.bam"


# Peak calling for bdgdiff pairwise comparisons

# Purpose:
        # This step runs MACS2 peak calling to generate narrowPeak and bedGraph files required for downstream pairwise comparisons using bdgdiff

#AF-RA

macs2 callpeak -B -f BAMPE -t $vl0029_25_AF_RA_H3K27ac -c $vl0037_25_AF_RA_Input -n 25_AF_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons 
macs2 callpeak -B -f BAMPE -t $vl0214_2231_AF_RA_H3K27ac -c $vl0339_2231_AF_RA_Input -n 2231_AF_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons 
macs2 callpeak -B -f BAMPE -t $vl0258_31_AF_RA_H3K27ac -c $vl0264_31_AF_RA_Input -n 31_AF_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons 
macs2 callpeak -B -f BAMPE -t $vl0260_313_AF_RA_H3K27ac -c $vl0266_313_AF_RA_Input -n 313_AF_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons 
macs2 callpeak -B -f BAMPE -t $vl0326_18_AF_RA_H3K27ac -c $vl0333_18_AF_RA_Input -n 18_AF_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons 
macs2 callpeak -B -f BAMPE -t $vl0033_28_AF_RA_H3K27ac -c $vl0041_28_AF_RA_Input -n 28_AF_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons

#SR-RA
macs2 callpeak -B -f BAMPE -t $vl0027_12b_SR_RA_H3K27ac -c $vl0035_12b_SR_RA_Input -n 12b_SR_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0031_27_SR_RA_H3K27ac -c $vl0039_27_SR_RA_Input -n 27_SR_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0218_33_SR_RA_H3K27ac -c $vl0234_33_SR_RA_Input -n 33_SR_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0261_20_SR_RA_H3K27ac -c $vl0267_20_SR_RA_Input -n 20_SR_RA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons

#AF-LA
macs2 callpeak -B -f BAMPE -t $vl0028_25_AF_LA_H3K27ac -c $vl0036_25_AF_LA_Input -n 25_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0032_28_AF_LA_H3K27ac -c $vl0040_28_AF_LA_Input -n 28_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0212_47_AF_LA_H3K27ac -c $vl0228_47_AF_LA_Input -n 47_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0215_31_AF_LA_H3K27ac -c $vl0231_31_AF_LA_Input -n 31_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0256_2231_AF_LA_H3K27ac -c $vl0262_2231_AF_LA_Input -n 2231_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0257_313_AF_LA_H3K27ac -c $vl0263_313_AF_LA_Input -n 313_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0327_18_AF_LA_H3K27ac -c $vl0334_18_AF_LA_Input -n 18_AF_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons

#SR-LA
macs2 callpeak -B -f BAMPE -t $vl0026_12b_SR_LA_H3K27ac -c $vl0034_12b_SR_LA_Input -n 12b_SR_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0213_45_SR_LA_H3K27ac -c $vl0229_45_SR_LA_Input -n 45_SR_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0217_20_SR_LA_H3K27ac -c $vl0340_20_SR_LA_Input -n 20_SR_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons
macs2 callpeak -B -f BAMPE -t $vl0328_27_SR_LA_H3K27ac -c $vl0335_27_SR_LA_Input -n 27_SR_LA --nomodel --extsize 301 -g 3.10e+09 --keep-dup all --outdir peak_calling_for_pairwise_comparisons




# Generate FE bdg files for visualization

#AF-RA
macs2 bdgcmp -t 25_AF_RA_treat_pileup.bdg -c 25_AF_RA_control_lambda.bdg -o 25_AF_RA_FE.bdg -m FE
macs2 bdgcmp -t 2231_AF_RA_treat_pileup.bdg -c 2231_AF_RA_control_lambda.bdg -o 2231_AF_RA_FE.bdg -m FE
macs2 bdgcmp -t 31_AF_RA_treat_pileup.bdg -c 31_AF_RA_control_lambda.bdg -o 31_AF_RA_FE.bdg -m FE
macs2 bdgcmp -t 313_AF_RA_treat_pileup.bdg -c 313_AF_RA_control_lambda.bdg -o 313_AF_RA_FE.bdg -m FE
macs2 bdgcmp -t 18_AF_RA_treat_pileup.bdg -c 18_AF_RA_control_lambda.bdg -o 18_AF_RA_FE.bdg -m FE
macs2 bdgcmp -t 28_AF_RA_treat_pileup.bdg -c 28_AF_RA_control_lambda.bdg -o 28_AF_RA_FE.bdg -m FE

#SR-RA
macs2 bdgcmp -t 12b_SR_RA_treat_pileup.bdg -c 12b_SR_RA_control_lambda.bdg -o 12b_SR_RA_FE.bdg -m FE
macs2 bdgcmp -t 27_SR_RA_treat_pileup.bdg -c 27_SR_RA_control_lambda.bdg -o 27_SR_RA_FE.bdg -m FE
macs2 bdgcmp -t 33_SR_RA_treat_pileup.bdg -c 33_SR_RA_control_lambda.bdg -o 33_SR_RA_FE.bdg -m FE
macs2 bdgcmp -t 20_SR_RA_treat_pileup.bdg -c 20_SR_RA_control_lambda.bdg -o 20_SR_RA_FE.bdg -m FE

#AF-LA
macs2 bdgcmp -t 25_AF_LA_treat_pileup.bdg -c 25_AF_LA_control_lambda.bdg -o 25_AF_LA_FE.bdg -m FE
macs2 bdgcmp -t 28_AF_LA_treat_pileup.bdg -c 28_AF_LA_control_lambda.bdg -o 28_AF_LA_FE.bdg -m FE
macs2 bdgcmp -t 47_AF_LA_treat_pileup.bdg -c 47_AF_LA_control_lambda.bdg -o 47_AF_LA_FE.bdg -m FE
macs2 bdgcmp -t 31_AF_LA_treat_pileup.bdg -c 31_AF_LA_control_lambda.bdg -o 31_AF_LA_FE.bdg -m FE
macs2 bdgcmp -t 2231_AF_LA_treat_pileup.bdg -c 2231_AF_LA_control_lambda.bdg -o 2231_AF_LA_FE.bdg -m FE
macs2 bdgcmp -t 313_AF_LA_treat_pileup.bdg -c 313_AF_LA_control_lambda.bdg -o 313_AF_LA_FE.bdg -m FE
macs2 bdgcmp -t 18_AF_LA_treat_pileup.bdg -c 18_AF_LA_control_lambda.bdg -o 18_AF_LA_FE.bdg -m FE

#SR-LA
macs2 bdgcmp -t 12b_SR_LA_treat_pileup.bdg -c 12b_SR_LA_control_lambda.bdg -o 12b_SR_LA_FE.bdg -m FE
macs2 bdgcmp -t 45_SR_LA_treat_pileup.bdg -c 45_SR_LA_control_lambda.bdg -o 45_SR_LA_FE.bdg -m FE
macs2 bdgcmp -t 20_SR_LA_treat_pileup.bdg -c 20_SR_LA_control_lambda.bdg -o 20_SR_LA_FE.bdg -m FE
macs2 bdgcmp -t 27_SR_LA_treat_pileup.bdg -c 27_SR_LA_control_lambda.bdg -o 27_SR_LA_FE.bdg -m FE
