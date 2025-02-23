################
# GWAS ANALYSIS 
################

# done by R.Verdugo (https://www.researchgate.net/profile/Ricardo-Verdugo-2), modified by J. Soto

# This script is to analyze GWAS from DGRP and link SNPs with genes using dgrp.fb557.annot

# This script need a variant annotation file for Drosophila melanogaster (FB5.57) http://dgrp2.gnets.ncsu.edu/data.html
# This script also need the file output "gwas.all.assoc" from the platform DGRP2 http://dgrp2.gnets.ncsu.edu/ 
# The files used in the DGRP2 platform include: 37_DGRP2.csv, 38_DGRP2.csv, 39_DGRP2.csv, 40_DGRP2.csv, ctmax_DGRP2.csv and z_DGRP2.csv

# example for producing miami plot of ctmax vs z

# Run each gwas.all.assoc file in this script to annotate the variants

# Read data
assoc <- read.table("gwas.all.assoc", header=T)

#to transform dgrp.fb549.annot.txt that contains all SNPs from FlyBase5.57 to object snptable
snptable <- read.csv("~/dgrp.fb557.annot.txt", header = FALSE, sep = "\t")

#to assign column names
colnames(snptable) <- c("ID", "ref", "SiteClass", "NN")

#to merge the SNPs of data (referred as assoc) with the database from FlyBase "snptable" by the SNP ID
dgrp.fb557.annot_gwas_all_assoc <- merge(assoc,snptable, by="ID")

##AVERAGE
#to select a subset of SNPs that are having a AvgMixedPval less or equal than 1e-05
SNPs_10e5_assoc <- subset(dgrp.fb557.annot_gwas_all_assoc, AvgMixedPval <= 1.00e-05)

#extract everything after SiteClass that is inside the brackets
genes<- gsub("^SiteClass\\[([^[]+)\\].*$", "\\1", SNPs_10e5_assoc$SiteClass)

#to split when there are more than one gene (FBgn) separated by ;
sites.ls<- strsplit(genes ,";")

#to avoid everything that is after the "|" symbol
genes.ls <- lapply(sites.ls, gsub, pattern="\\|.*$", replacement="")

#to extract name of gene
name.ls<- lapply(sites.ls, gsub, pattern="^\\w*\\|(.*)\\|\\w*\\|.*$", replacement="\\1")

#to extract site of SNP
place.ls <- lapply(sites.ls, gsub, pattern=".*\\|(\\w*)\\|.*$", replacement="\\1")

#to keep the SNP name
names(genes.ls) <- csv$ID
genes.ids <- sub("\\d$", "", names(unlist(genes.ls)))

#to merge the genes list ID with the FBgns
rs_gene <- cbind(RS=genes.ids, Gene= unlist(genes.ls), Name=unlist(name.ls), Site=unlist(place.ls))


# This produce a data.frame with the annotated SNPs significantly associated with the phenotype
# to save it as table separated by ",", CTmax used as example:
write.table(rs_gene, file="57_Heat_tolerance_CTmax_SNPs_10e5_rsbygene.csv", row.name=F, col.names=T, sep=",", quote=F)

#############
# MIAMI PLOTS
#############
# this plot shows 2 Manhattan plots agains each other. Ctmax and Z used as examples

library(qqman)

# A CTmax
# Load the "gwas.all.assoc" output obtained from the platform DGRP2, using the CTmax data
assocA <- read.table("gwas.CTmax.all.assoc", header=T)

assocA$CHR <- sub ("_.*_.*$", "", assocA$ID)
assocA$BP <- sub ("^.*_(.*)_.*$", "\\1", assocA$ID)
assocA$CHR <- as.numeric(as.factor(assocA$CHR))
assocA$BP <- as.numeric(assocA$BP)

Mixed_GxEA <- cbind(assocA[,c("CHR", "BP")], SNP=assocA$ID, P=assocA$AvgMixedPval)

# Load table obtained before with the SNPs significantly associated with the phenotype, for CTmax
snpA <- read.table("57_Heat_tolerance_CTmax_SNPs_10e5_rsbygene.csv", sep = ",", header = T)
#view max p in plot
summary(-log10(Mixed_GxEA$P))

# B Z
# Load the "gwas.all.assoc" output obtained from the platform DGRP2, using the Z data
assocA <- read.table("gwas.z.all.assoc", header=T)

assocB$CHR <- sub ("_.*_.*$", "", assocB$ID)
assocB$BP <- sub ("^.*_(.*)_.*$", "\\1", assocB$ID)
assocB$CHR <- as.numeric(as.factor(assocB$CHR))
assocB$BP <- as.numeric(assocB$BP)

Mixed_GxEB <- cbind(assocB[,c("CHR", "BP")], SNP=assocB$ID, P=assocB$MaleMixedPval)

# Load table obtained before with the SNPs significantly associated with the phenotype, for Z
snpB <- read.table("57_Heat_tolerance_z_SNPs_10e5_rsbygene.csv", sep = ",", header = T)
#view max p in plot
summary(-log10(Mixed_GxEB$P))

# PLOT
# Make the miami plot and save it in .png format
png(file="Miami_ctmax_z_avg.png",width=3000, height=2000, res = 300)
par(mfrow=c(2,1))
par(mar=c(1,5,1,3))
manhattan(na.omit(Mixed_GxEA),chrlabs = c("2L", "2R", "3L", "3R", "4", "X"), genomewideline=FALSE, highlight=as.vector(snpA[,1]),yaxt = "n",ylim=c(0,10))
axis(side = 2,las=2,at=c(0,1,2,3,4,5,6,7,8,9,10))
legend("topright",bty = "n", legend = "CTmax: 76 SNPs",cex = 2)
par(mar=c(1,5,1,3))
manhattan(na.omit(Mixed_GxEB),chrlabs = c("2L", "2R", "3L", "3R", "4", "X"), genomewideline=FALSE, highlight=as.vector(snpB[,1]),yaxt = "n",ylim=c(10,0),xlab="",xaxt="n")
axis(side = 2,las=2,at=c(0,1,2,3,4,5,6,7,8,9,10))
legend("bottomright",bty = "n", legend = "Z: 140 SNPs",cex = 2)
dev.off()