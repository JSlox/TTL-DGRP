####################
# EFFECT SIZE PLOT
####################

library(reshape2)
library(ggplot2)
library(ggplotify)
library(gridExtra)
library(ggpubr)

#CTmax --------------------------------------------------------------------------
CTmax_genot<-read.table("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/CTmax_snp_geno.csv", header = T,sep =";",check.names = F)
CTmax_snp_genes <- read.table("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/CTmax_SNP.csv", header = T,sep =";")

for (i in 3:ncol(CTmax_genot)){
  CTmax_genot[,i] <- as.factor(CTmax_genot[,i]) 
}

CTmax_genot_melt <- melt(CTmax_genot[,-1])
lista_graficos <- list()

for (i in 1:8){
  propiedades_snp <- subset(CTmax_snp_genes,ID==colnames(CTmax_genot_melt[i]))
  grafico <- ggplot(subset(CTmax_genot_melt,!is.na(CTmax_genot_melt[,i])), aes(x=subset(CTmax_genot_melt[,i],!is.na(CTmax_genot_melt[,i])), y= value, color=subset(CTmax_genot_melt[,i],!is.na(CTmax_genot_melt[,i])))) +
    geom_boxplot() + theme(panel.background = element_blank(),
                           legend.position = "none", 
                           panel.border = element_rect(colour = "black", fill=NA))+
    xlab(paste(colnames(CTmax_genot_melt[i]),propiedades_snp$Name,sep="\n")) +
    ylab("")+
    annotate("text", x=2.2, y=54, label= propiedades_snp$AvgEff) + 
    annotate("text", x =2.2, y=53, label = propiedades_snp$AvgMixedPval)+
    scale_color_brewer(palette="Dark2")+
    ylim(44,54)
  grob_grafico <- as.grob(grafico)
  nombre_grafico <- paste0("SNP",i)
  lista_graficos[[nombre_grafico]] <- grob_grafico
}

effect_ctmax_plot <- grid.arrange(grobs=c(c(lista_graficos[1:8])), ncol=4,left=("CTmax(°C)"))

#Z ------------------------------------------------------------------------------
z_genot <- read.table("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/Z_snp_geno.csv", header = T,sep =";",check.names = F)
z_snp_genes <- read.table("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/Z_SNP.csv", header = T,sep =";")

for (i in 3:ncol(z_genot)){
  z_genot[,i] <- as.factor(z_genot[,i]) 
}

z_genot_melt <- melt(z_genot[,-1])
lista_graficos <- list()

for (i in 1:8){
  propiedades_snp <- subset(z_snp_genes,ID==colnames(z_genot_melt[i]))
  grafico <- ggplot(subset(z_genot_melt,!is.na(z_genot_melt[,i])), aes(x=subset(z_genot_melt[,i],!is.na(z_genot_melt[,i])), y= value, color=subset(z_genot_melt[,i],!is.na(z_genot_melt[,i])))) +
    geom_boxplot() + theme(panel.background = element_blank(),
                           legend.position = "none", 
                           panel.border = element_rect(colour = "black", fill=NA))+
    xlab(paste(colnames(z_genot_melt[i]),propiedades_snp$Name,sep="\n")) +
    ylab("")+
    annotate("text", x=2.2, y=11, label= propiedades_snp$AvgEff) + 
    annotate("text", x =2.2, y=10, label = propiedades_snp$AvgMixedPval)+
    scale_color_brewer(palette="Dark2")+
    ylim(3,11)
  grafico <- grafico + scale_y_continuous(breaks=c(4,6,8,10), labels=c("4","6","8","10"))
  grob_grafico <- as.grob(grafico)
  nombre_grafico <- paste0("SNP",i)
  lista_graficos[[nombre_grafico]] <- grob_grafico
}

effect_z_plot <- grid.arrange(grobs=c(c(lista_graficos[1:8])), ncol=4,left=("z(°C)"))


png(file="Effect_box.png",width=3000, height=2500, res = 300)
ggarrange(effect_ctmax_plot, effect_z_plot, labels=c("A","B"), ncol=1, nrow=2)
dev.off()