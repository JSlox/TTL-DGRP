####################
# EFFECT SIZE PLOT
####################

#LOAD REQUIRED PACKAGES
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(reshape2, ggplot2, ggplotify, gridExtra, ggpubr)

#CTmax --------------------------------------------------------------------------
# Load data of CTmax and genotype per snp and snps effect size information 
CTmax_genot<-read.table("CTmax_snp_geno.csv", header = T,sep =";",check.names = F)
CTmax_snp_genes <- read.table("CTmax_SNP.csv", header = T,sep =";")

# prepare data
for (i in 3:ncol(CTmax_genot)){
  CTmax_genot[,i] <- as.factor(CTmax_genot[,i]) 
}

CTmax_genot_melt <- melt(CTmax_genot[,-1])
# make empty list to store graphs
lista_graficos <- list()

# make box plot per snp and store it in list
for (i in 1:8){
  propiedades_snp <- subset(CTmax_snp_genes,ID==colnames(CTmax_genot_melt[i]))
  grafico <- ggplot(subset(CTmax_genot_melt,!is.na(CTmax_genot_melt[,i])), aes(x=subset(CTmax_genot_melt[,i],!is.na(CTmax_genot_melt[,i])), y= value, color=subset(CTmax_genot_melt[,i],!is.na(CTmax_genot_melt[,i])))) +
    geom_boxplot() + theme(panel.background = element_blank(),
                           axis.text.y= element_text(colour="black"), axis.ticks = element_line(color = "black"),
                           axis.text.x= element_text(colour="black"), legend.position = "none", 
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

# arrange graphs
effect_ctmax_plot <- grid.arrange(grobs=c(c(lista_graficos[1:8])), ncol=4,left=("CTmax(°C)"))

#Z ------------------------------------------------------------------------------
z_genot <- read.table("Z_snp_geno.csv", header = T,sep =";",check.names = F)
z_snp_genes <- read.table("Z_SNP.csv", header = T,sep =";")

for (i in 3:ncol(z_genot)){
  z_genot[,i] <- as.factor(z_genot[,i]) 
}

z_genot_melt <- melt(z_genot[,-1])
lista_graficos <- list()

for (i in 1:8){
  propiedades_snp <- subset(z_snp_genes,ID==colnames(z_genot_melt[i]))
  grafico <- ggplot(subset(z_genot_melt,!is.na(z_genot_melt[,i])), aes(x=subset(z_genot_melt[,i],!is.na(z_genot_melt[,i])), y= value, color=subset(z_genot_melt[,i],!is.na(z_genot_melt[,i])))) +
    geom_boxplot() + theme(panel.background = element_blank(),axis.text.y= element_text(colour="black"), axis.ticks = element_line(color = "black"),
                           axis.text.x= element_text(colour="black"),
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

# arrange graphs for CTmax and z
effect_ctmax_z <- ggarrange(effect_ctmax_plot, effect_z_plot, labels=c("A","B"), ncol=1, nrow=2)
ggsave(file="effect_ctmax_z.svg", plot=effect_ctmax_z, width=3000, height=2500, dpi = 300, units="px")
