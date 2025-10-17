

##########################################
##########################################
###### HOMER TF Enrichment Analysis ######
##########################################
##########################################

## The findMotifsGenome.pl script (-size given) from the HOMER software package (v5.1) was used to test for enrichment of TF motifs at H3K27ac-enriched regions and DMRs. 
## For H3K27ac data, the sets of cellective peaks from each category were used as background to identify disease or anatomical side-specific signal. 
## For DMRs, all candidate regions with aggregated CpG methylation signal were used as background. 

cd Rodriguez_2025_AFepigenome/TF_enrichment_analysis

# Install HOMER
perl configureHomer.pl -install homer
# Add HOMER programs to executable path

# Install hg38
perl configureHomer.pl -install hg38


# Run HOMER from bin folder

## Files with set of regions need to be edited to make them compatible with HOMER input:
    # Column1: chromosome (with "chr" prefix)
    # Column2: starting position
    # Column3: ending position
    # Column5: Strand ()
    # E.g.: 
    # sed 's/$/\t /' < AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr.bed > AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_edited1.bed
    # sed 's/$/\t+/' < AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_edited1.bed > AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_edited2.bed


#### ChIP-seq H3K37ac data

# AF-LA
findMotifsGenome.pl AF_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_c2_edited2.bed hg38 AFLA_MotifOutput_GenrichBackground/ -size given -mask -bg AF_LA_Genrich_greylistExcluded_g10-edited3.bed

# AF-RA
findMotifsGenome.pl AF_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_c2_edited2.bed hg38 AFRA_MotifOutput_GenrichBackground/ -size given -mask -bg AF_RA_Genrich_greylistExcluded_g10-edited3.bed

# SR-LA
findMotifsGenome.pl SR_LA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_c2_edited2.bed hg38 SRLA_MotifOutput_GenrichBackground/ -size given -mask -bg SR_LA_Genrich_greylistExcluded_g10-edited3.bed

# SR-RA
findMotifsGenome.pl SR_RA_enriched_nonoverlappingwithallsets_GenrichGreylistedFiltered_chr_c2_edited2.bed hg38 SRRA_MotifOutput_GenrichBackground/ -size given -mask -bg SR_RA_Genrich_greylistExcluded_g10-edited3.bed


#### DNAm data - Hypomethylated and hypermethylated DMRs

# AF-hypermethylated
findMotifsGenome.pl dmrs_dmrff_betanorm_Chris_DMRs_hypermeth_editedforHOMER.bed hg38 DMRs_hypermeth/ -size given -mask -bg dmrs_dmrff_betanorm_Chris_DMRs_HOMER.bed

# AF-hypomethylated
findMotifsGenome.pl dmrs_dmrff_betanorm_Chris_DMRs_hypometh_editedforHOMER.bed hg38 DMRs_hypometh/ -size given -mask -bg dmrs_dmrff_betanorm_Chris_DMRs_HOMER.bed


