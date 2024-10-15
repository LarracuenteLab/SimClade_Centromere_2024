# SimClade_Centromere_2024
This repository contains code required to analyze data in Courret et al. 2024 (https://doi.org/10.1101/2023.08.22.554357). In this manuscript, we study centromere evolution in the Drosophila simulans clade – D. simulans, D. sechellia, and D. mauritiana – compared to their relative, D. melanogaster. We use chromatin profiling (CUT&Tag) combined with high resolution fluorescence in situ hybridization on stretched DNA, and evolutionary analyses to characterize and compare all centromeres across these species.

The directories contain the R code and the intermediate files necessary to make the Figures, and the bash script used to analyze the raw data. R scripts are organized by Figure type. The subdirectories contain code for three types of plots and the data analysis pipeline as follows:

CoveragePlot

Description: This directory contains all the R scripts necessary to make Figures 1ADG-2-3-4AB-5-S1-S2-S2-S5. The R script require the library karyoploteR, regioneR, GenomicRanges, rtracklayer, IRanges, devtools, stringr. This folder contain a subfolder named "Reference" that contain the genome object (chromosome.size and "cytoband") necessary to make those figures. 

IDR

Description: This directory contains the R scripts necessary to make Figure S4. This directory also contains the file necessary to make this figure: the chrome.size object and the q0.Top-idr (tab delimited file containing the top peaks conserved between replicates).

TE_enrichment_Plot:

Description: This directory contains the R script necessary to make Figure 1BEH. This directory also contains the file necessary to make this figure: the "CID_RPM_Top_newname.out" is tab delimited file that contain the top most enriched TE.

Pipeline 

Description : This directory contains the bash script used to analyzed the CUT&Tag data.

script: CutandTag.sh - performs the analysis of the CUT&Tag raw data to obtain the coverage plots. (for Figure 1ADG-2-3-4AB-5-S1-S2-S2-S5)

script: CutandTag_TEcount.sh - estimates the TE enrichment (for Figure 1BEH)

script: CutandTag_peakCalling.sh - does the peak calling (for Figure 1ADG-2-3-4AB-5-S1-S2-S2-S5).

script: Blast_Sat.sh - identifies the satellite locations in the genome using Blast (for Figure S7).
