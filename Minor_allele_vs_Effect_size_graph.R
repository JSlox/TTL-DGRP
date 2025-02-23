#############################################
# Plot Minor allele effect vs effect of SNP
#############################################

library(scales)

# LOAD AND PREPARE DATA ------------------------------------------------------
#CTmax
## load results from gwas.top.annot file
effect_CTmax <- read.table("ctmax_gwas.top.annot", sep = " ", header = T)

## subset data (AvgMixedPval, DiffMixedPval, FemaleMixedPval or MaleMixedPval)
sub_effect_CTmax <- subset(effect_CTmax, AvgMixedPval <= 1.00e-05)
sub_effect_CTmax <- sub_effect_CTmax[,c("ID","MAF","AvgEff")]
#write.table(sub_effect_CTmax,"CTmax_effect.csv")

## same but for z
effect_z <- read.table("z_gwas.top.annot", sep = " ", header = T)
sub_effect_z <- subset(effect_z, AvgMixedPval <= 1.00e-05)
sub_effect_z <- sub_effect_z[,c("ID","MAF","AvgEff")]
#write.table(sub_effect_z,"Z_effect.csv")

# PLOT ----------------------------------------------------------------------
## Plot frequency of minor allele vs effect of SNP
## Change plot size and axis X acording to range of effect
## Change color of points too
## Change name of file
par(mar=c(1,1,1,1))
tiff(file="Avg_Effect_CTmax_z.tiff",units="px",width=2100, height=2000, res=300)
layout(mat = matrix(c(1, 2), nrow = 2, ncol = 1), # matriz and order of plots
       heights = c(3.5, 1), # Heights of the rows
       widths = c(1, 1)) # Widths of the columns
#PLOT 1 EFFECT
par(mar=c(3.5,4,1,1))
plot(c(-1,1), c(0, 0.5), type = "n", axes = FALSE, xlab = "", ylab = "")
for (i in 1:nrow(sub_effect_CTmax)) {
points(sub_effect_CTmax[i,3], sub_effect_CTmax[i,2], type = "p", cex = 2, lwd=2 , pch = 21, col = alpha("black", 0.7), bg=alpha("green3", 0.5))
}

for (i in 1:nrow(sub_effect_z)) {
  points(sub_effect_z[i,3], sub_effect_z[i,2], type = "p", cex = 2, lwd=2 , alpha=0.5, pch = 21, col = alpha("black", 0.7), bg=alpha("grey", 0.5))
}

axis(side = 1, at = c(-1.0,-0.5,0.0,0.5,1.0), labels = c("-1.0","-0.5","0.0","0.5","1.0"), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(0.1,0.2,0.3,0.4,0.5), labels = c("0.1","0.2","0.3","0.4","0.5"), cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(xlab = "Effect size", cex.lab = 2, mgp = c(2.5, 0, 0))
title(ylab = "Minor allele frequency", cex.lab = 2, mgp = c(2.5, 0, 0))
abline(v=0, col="black")
legend("bottomright", bty = "n", pch = 21, col = c("black","black"), pt.cex = 1.5, pt.bg = c("green3","grey"), legend = c("CTmax", "z"), text.col = c("green3","grey"), x.intersp = 0.5, y.intersp = 1, cex = 1.7)

# PLOT 2 ARROWS
par(mar=c(0.5,4,0.5,1))#The default is c(5.1, 4.1, 4.1, 2.1) bottom, left, top, and right
plot(c(-1,1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
arrows(x0=-1, y0=0.9, x1=1, y1=0.9, length = 0.25, angle = 30,
       code = 3, col = "green3", lty = "solid", lwd = 8) # lty=type line, lwd=width line
arrows(x0=-1, y0=0.4, x1=1, y1=0.4, length = 0.25, angle = 30,
       code = 3, col = "grey", lty = "solid", lwd = 8) # lty=type line, lwd=width line
text(x=-0.9,y=0.65,labels=paste("Minor allele","increases CTmax",sep="\n"),pos=4,cex=1.2)
text(x=0.9,y=0.65,labels=paste("Minor allele","decreases CTmax",sep="\n"),pos=2,cex=1.2)
text(x=-0.9,y=0.15,labels=paste("Minor allele","decreases sensitivity",sep="\n"),pos=4,cex=1.2)
text(x=0.9,y=0.15,labels=paste("Minor allele","increases sensitivity",sep="\n"),pos=2,cex=1.2)
dev.off()
