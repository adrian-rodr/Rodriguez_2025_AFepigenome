# Rodriguez_2025_AFepigenome

This archive contains processed files and scripts associated with **Rodriguez et al (2025)**, *"An epigenomic investigation of atrial fibrillation in a matched left and right atrial human cohort"*. doi: https://doi.org/10.1101/2025.08.29.673028

The repository includes data derived from genome-wide sequencing and microarray analysis, as well as code used for downstream analyses and for generating the figures presented in the publication.

The raw ChIP-sequencing and DNA methylation array data reported in this study have been deposited in the GEO under accession numbers GSE227793 and GSE291249, respectively.

## ChIP-seq data analysis

*Basic processing: alignment and peak calling*

The script [`Read_alignment_AF-ChIP-seq_dataset.pl`] was used to submit alignment jobs. [`BAM_files_preprocessing.pl`]performs basic preprocessing of BAM files prior to peak calling (i.e., removing multi-mapping reads and subsampling to 20M reads). Peak calling was conducted using [`MACS2_peak_calling.sh`], which generate **.bedGraph** and **narrowPeak** files used for downstream analyses. This script also produces fold-enrichment (FE) **.bdg** files, which can be loaded into genome browsers for visualisation of peak regions.

**Note**
The above scripts are provided as-is for transparency. They were written for execution on the QMUL Apocrita HPC cluster and may not reproduce identically in other computing environments. Users intending to adapt these scripts should update cluster-specific resources (e.g. genome indices, module load path, etc.). 

*Definition of regulatory regions across each sample type*

## Definition of enriched H3K27ac regions across sample groups

## RNA-seq analysis of differential gene expression

## RT-qPCR and ChIP-qPCR

## DNAm EPICv2 array analysis

## Transcription Factor Enrichment analysis