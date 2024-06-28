# TTL-DGRP

Scripts to analyze thermal tolerance landscape of 100 lines of the DGRP as shown in Soto et al. 2024.

A summary of the various inputs and script needed:   
INPUT ----- > SCRIPT   
KO_4temp_100DGRP.xlsx -----> Reaction_norms.R.; TDT_PLOT.R; LMMs.R   
ctmax_gwas.top.annot; z_gwas.top.annot -----> Minor_allele_vs_Effect_size_graph.R   
CTmax_snp_geno.csv; CTmax_SNP.csv; Z_snp_geno.csv; Z_SNP.csv -----> Effect_size_plot.R   
RNAi_validation.xlsx -----> RNAi_validation.R   
Genetic_Correlation.txt -----> Plot_genetic_correlation.R   
37_DGRP2.csv; 38_DGRP2.csv; 39_DGRP2.csv; 40_DGRP2.csv; ctmax_DGRP2.csv; z_DGRP2.csv -----> GWAS_analysis.R   

## Reaction norms and TDT analysis
Data file ***KO_4temp_100DGRP.xlsx*** contains all data from the thermal assays performed on 100 lines of the DGRP over 4 temperatures, indicating cell of the fly (celda), experimental block (bloque), date (fecha), sex, genotype ID (ID), code of genotype (Geno), experimental temperature (temp) and knockdown time on seconds (KO). This is the input file for 3 scripts:

  - (1) ***Reaction_norms.R***, that makes the reaction norms figures (1 and 2).
  - (2) ***TDT_PLOT.R***, that make the plots of the TDT lines (Figure 3).
  - (3) ***LMMs.R***, that performs the linear mixed models on Knockdown data and CTmax and z data (Tables 1,2 and 3, Table S1, S3 and S5).

Script ***Thermal_landscape.R*** contains a function to calculate thermal tolerance landscape parameters (CTmax and Z).

## GWAS analysis
The ***GWAS_analysis.R*** script annotates genomic variants output of the DGRP2 platform and also performs miami plots for the GWAS analysis. This script needs various inputs that includes the output of the DGRP2 platform and a Variant annotation file for drosophila melanogaster. The files ***37_DGRP2.csv, 38_DGRP2.csv, 39_DGRP2.csv, 40_DGRP2.csv, ctmax_DGRP2.csv and z_DGRP2.csv*** all are analyzed separately in the DGRP2 platform, each containing the knockdown time mean per line and sex for each temperature (37-40), or the ctmax and z mean per line and sex for the ctmax and z files. Further instructions in script (Figure 4 and S2-S7, Tables S6-S11)

## Effect size vs minor allele frequency
Files ***ctmax_gwas.top.annot*** and ***z_gwas.top.annot*** contains the SNPs significantly associated with CTmax and z, respectively. This files are the input on the ***Minor_allele_vs_Effect_size_graph.R***, that plot the minor allele frequency vs the effect size for each of the SNPs (Figure 5).

## Alleles effect
Files ***CTmax_snp_geno.csv*** and ***CTmax_SNP.csv*** contains information of 8 of the SNPs significantly associated with CTmax. ***CTmax_snp_geno.csv*** indicates the CTmax phenotypic value for each line of the DGRP and the genotype of the line for each SNP (mayor allele = 0 and minor allele = 2). ***CTmax_SNP.csv*** indicates the ID if the SNPs, gene associated, effect size and significance of the association with CTmax. Files ***Z_snp_geno.csv*** and ***Z_SNP.csv*** contains the same but for z. All of this files are the input for the ***Effect_size_plot.R*** (Figure 6).

## RNAi validation
File ***RNAi_validation.xlsx*** contains the results of the experiments to validate the phenotypic effect of some of the candidate genes (mam, shot, robo3, KCNQ) using RNAi. indicating cell of the fly (celda), experimental replicate (Replica), date (fecha), sex, stock ID of fly (ID), code of genotype (Geno), experimental temperature (temp) and knockdown time on seconds (KO). This is the input for the ***RNAi_validation.R*** script, that performs the statistical analysis of the results and the figures comparing treatment and control (Figure 7 and table S14).

## Genetic correlation
File ***Genetic_Correlation.txt*** contains the input for the MTDFREML software and the ***Plot_genetic_correlation.R*** script, indicating knockdown time for each fly in each temperature. The script ***Plot_genetic_correlation.R*** plots the genetic correlation across temperatures (Figure S1 and Table 3).