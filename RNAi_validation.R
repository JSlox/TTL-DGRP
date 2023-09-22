library(lme4) # linear mixed model
library(lmerTest) # singnificance of random effects
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(rstatix)
library(dplyr)
library(gridExtra)
library(ggplotify)

################
# RNAi ANALYSIS 
################

# READ DATA -------------------------------------------------------------------

library(readxl)
library(reshape2)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tidyverse)
library(gridExtra)
setwd("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data")
data <- read_excel("RNAi_validation.xlsx")

data$Geno <- as.factor(data$Geno)
data$temp <- as.factor(data$temp)
data$sex <- as.factor(data$sex)
data$Replica <- as.factor(data$Replica)
data$Fecha <- as.factor(data$Fecha)
data$Celda <- as.factor(data$Celda)
data$KO <- as.numeric(data$KO/60)

datarobo3 <- droplevels(subset(data,Geno=="robo3-attp40"|Geno=="CONTROL-attp40"))
datamam <- droplevels(subset(data,Geno=="mam-kk"|Geno=="CONTROL-kk"))
datashot <- droplevels(subset(data,Geno=="shot-sh"|Geno=="CONTROL-sh"))
dataKCNQ2 <- droplevels(subset(data,Geno=="KCNQ-attp2"|Geno=="CONTROL-attp2"))


# THERMAL LANDSCAPE --------------------------------

# CTMAX AND Z
source("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/TTL-DGRP/Thermal_landscape.R")

R1_data <- droplevels(subset(data,Replica=="R1"))
TDT_table_R1 <- TDT(R1_data)
Replica <- rep("R1",16)
TDT_table_R1 <- cbind(TDT_table_R1,Replica)

R2_data <- droplevels(subset(data,Replica=="R2"))
TDT_table_R2 <- TDT(R2_data)
Replica <- rep("R2",16)
TDT_table_R2 <- cbind(TDT_table_R2,Replica)

R3_data <- droplevels(subset(data,Replica=="R3"))
TDT_table_R3 <- TDT(R3_data)
Replica <- rep("R3",16)
TDT_table_R3 <- cbind(TDT_table_R3,Replica)

TDT_full <- rbind(TDT_table_R1,TDT_table_R2,TDT_table_R3)
TDT_full$sex <- as.factor(TDT_full$sex)

##sexes comparison of CTmax and z

#CTmax
shapiro.test(TDT_full$CTmax)
fligner.test(TDT_full$CTmax,TDT_full$sex)
wilcox.test(subset(TDT_full, sex=="F")$CTmax, subset(TDT_full, sex=="M")$CTmax, alternative = "two.sided")

#Z
shapiro.test(TDT_full$Z)
fligner.test(TDT_full$Z,TDT_full$sex)
wilcox.test(subset(TDT_full, sex=="F")$Z, subset(TDT_full, sex=="M")$Z, alternative = "two.sided")

cor(TDT_full$CTmax,TDT_full$Z)
anova(manova(cbind(CTmax,Z)~sex, data= TDT_full))#no se cumplen supuestos
summary(manova(cbind(CTmax,Z)~sex, data= TDT_full))# no se cumplen supuestos

# Boxplot CTmax
png(file="Box_CTmax.png",width=3000, height=1700, res = 300)
plot(c(1, 4.05), c(0, 200), type = "n", axes = FALSE, xlab = "", ylab = "")
boxplot(TDT_table$CTmax~TDT_table$sex, xlab = "", ylab="", 
        col = brewer.pal(9, "Set1")[c(3, 1, 2,9)],  cex.axis=2, cex.lab=2.2, bty="n" )
title(ylab = "CTmax", cex.lab = 2.2, mgp = c(2.6, 0, 0))
title(xlab = "Sexo", cex.lab = 2.2, mgp = c(3.5, 0, 0))
dev.off()

# Boxplot Z
png(file="Box_z.png",width=3000, height=1700, res = 300)
plot(c(1, 4.05), c(0, 200), type = "n", axes = FALSE, xlab = "", ylab = "")
boxplot(TDT_table$Z~TDT_table$sex, xlab = "", ylab="", 
        col = brewer.pal(9, "Set1")[c(3, 1, 2,9)],  cex.axis=2, cex.lab=2.2, bty="n" )
title(ylab = "Z", cex.lab = 2.2, mgp = c(2.6, 0, 0))
title(xlab = "Sexo", cex.lab = 2.2, mgp = c(3.5, 0, 0))
dev.off()

# TDT PLOT
rep.statics <- 10^with(data,tapply(log10(KO),list(temp,Geno,sex),mean))	#geometric mean for each replicate
out <- matrix(,(nlevels(data$Geno)*8),1)
for(i in 1:(nlevels(data$Geno)*8)){out[i,1] <- rep.statics[(i)]}
colnames(out) <- c("TimeKO")
sex <- factor(rep(c(1,2),each=(nlevels(data$Geno)*4)),labels=c("females","males"))
Geno <- c(levels(data$Geno))
Geno <- factor(rep(c(1:(nlevels(data$Geno))),each=4,2),labels=Geno)
temp <- factor(rep(c(37:40),(nlevels(data$Geno)*2)))
statics.rep <- data.frame(temp,Geno,sex,out)


### TDT curve plot
rep.statics <- 10^with(data,tapply(log10(KO),list(temp,Geno,sex),mean))	#geometric mean for each replicate
data.F <- statics.rep[which(statics.rep[,"sex"]=="females"),]
data.M <- statics.rep[which(statics.rep[,"sex"]=="males"),]
head(data.F)
head(data.M)

# Females
png(file="TDThembras.png",width=3000, height=1700, res = 300)
plot(c(1, 4.05), c(0, 120),xlab="",ylab="",axes = FALSE,log="y",las=1,ylim=c(1,120),xlim=c(36,46),xaxs="i",yaxs="i")
for(i in 1:9){
  abline(lm(log10(data.F[(i*4-3):(i*4),4]) ~ c(37:40)),col="firebrick1", lwd=1)} 

axis(side = 1, at = c(36,38,40,42,44,46,48,50,52,54), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(1,5,20,100), cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Tiempo de colapso (min)", cex.lab = 2.2, mgp = c(2.5, 0, 0))
title(xlab = "Temperatura (°C)", cex.lab = 2.2, mgp = c(3.5, 0, 0))
dev.off()

# Males
png(file="TDTmachos.png",width=3000, height=1700, res = 300)
plot(c(1, 4.05), c(0, 120),xlab="",ylab="",axes = FALSE,log="y",las=1,ylim=c(1,120),xlim=c(36,48),xaxs="i",yaxs="i")
for(i in 1:9){
  abline(lm(log10(data.M[(i*4-3):(i*4),4]) ~ c(37:40)),col="dodgerblue1", lwd=1)} 

axis(side = 1, at = c(36,38,40,42,44,46,48,50,52,54), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(1,5,20,100), cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Tiempo de colapso (min)", cex.lab = 2.2, mgp = c(2.5, 0, 0))
title(xlab = "Temperatura (°C)", cex.lab = 2.2, mgp = c(3.5, 0, 0))
dev.off()


# LMM ON CTMAX AND Z -----------------------------------

#CTmax y Z

TDT_datarobo3 <- droplevels(subset(TDT_full,Geno=="robo3-attp40"|Geno=="CONTROL-attp40"))
TDT_datamam <- droplevels(subset(TDT_full,Geno=="mam-kk"|Geno=="CONTROL-kk"))
TDT_datashot <- droplevels(subset(TDT_full,Geno=="shot-sh"|Geno=="CONTROL-sh"))
TDT_dataKCNQ2 <- droplevels(subset(TDT_full,Geno=="KCNQ-attp2"|Geno=="CONTROL-attp2"))
TDT_dataKCNQ40 <- droplevels(subset(TDT_full,Geno=="KCNQ-attp40"|Geno=="CONTROL-attp40"))

# CTmax
CTmax.robo3 <- lm(CTmax ~ Replica+Geno, data = subset(TDT_datarobo3,sex=="M"))
summary(CTmax.robo3)
anova(CTmax.robo3)
wilcox.test (subset(TDT_datarobo3,sex=="M")$CTmax~subset(TDT_datarobo3,sex=="M")$Geno,alternative = "two.sided", paired=F)

CTmax.mam <- lm(CTmax ~ Replica+Geno, data = subset(TDT_datamam,sex=="M"))
summary(CTmax.mam)
anova(CTmax.mam)
wilcox.test (subset(TDT_datamam,sex=="M")$CTmax~subset(TDT_datamam,sex=="M")$Geno,alternative = "two.sided", paired=F)

CTmax.shot <- lm(CTmax ~ Replica+Geno, data = subset(TDT_datashot,sex=="M"))
summary(CTmax.shot)
anova(CTmax.shot)
wilcox.test (subset(TDT_datashot,sex=="F")$CTmax~subset(TDT_datashot,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_datashot,sex=="M")$CTmax~subset(TDT_datashot,sex=="M")$Geno,alternative = "two.sided", paired=F)

CTmax.KCNQ2 <- lm(CTmax ~ Replica+Geno, data = subset(TDT_dataKCNQ2,sex=="F"))
summary(CTmax.KCNQ2)
anova(CTmax.KCNQ2)
wilcox.test (subset(TDT_dataKCNQ2,sex=="F")$CTmax~subset(TDT_dataKCNQ2,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_dataKCNQ2,sex=="M")$CTmax~subset(TDT_dataKCNQ2,sex=="M")$Geno,alternative = "two.sided", paired=F)

CTmax.KCNQ40 <- lm(CTmax ~ Replica+Geno, data = subset(TDT_dataKCNQ40,sex=="M"))
summary(CTmax.KCNQ40)
anova(CTmax.KCNQ40)
wilcox.test (subset(TDT_dataKCNQ40,sex=="F")$CTmax~subset(TDT_dataKCNQ40,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_dataKCNQ40,sex=="M")$CTmax~subset(TDT_dataKCNQ40,sex=="M")$Geno,alternative = "two.sided", paired=F)

# Z
Z.robo3 <- lm(Z ~ Replica+Geno, data = subset(TDT_datarobo3,sex=="F"))
summary(Z.robo3)
anova(Z.robo3)
wilcox.test (subset(TDT_datarobo3,sex=="F")$Z~subset(TDT_datarobo3,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_datarobo3,sex=="M")$Z~subset(TDT_datarobo3,sex=="M")$Geno,alternative = "two.sided", paired=F)

Z.mam <- lm(Z ~ Replica+Geno, data = subset(TDT_datamam,sex=="M"))
summary(Z.mam)
anova(Z.mam)
wilcox.test (subset(TDT_datamam,sex=="F")$Z~subset(TDT_datamam,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_datamam,sex=="M")$Z~subset(TDT_datamam,sex=="M")$Geno,alternative = "two.sided", paired=F)

Z.shot <- lm(Z ~ Replica+Geno, data = subset(TDT_datashot,sex=="F"))
summary(Z.shot)
anova(Z.shot)
wilcox.test (subset(TDT_datashot,sex=="F")$Z~subset(TDT_datashot,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_datashot,sex=="M")$Z~subset(TDT_datashot,sex=="M")$Geno,alternative = "two.sided", paired=F)

Z.KCNQ2 <- lm(Z ~ Replica+Geno, data = subset(TDT_dataKCNQ2,sex=="M"))
summary(Z.KCNQ2)
anova(Z.KCNQ2)
wilcox.test (subset(TDT_dataKCNQ2,sex=="F")$Z~subset(TDT_dataKCNQ2,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_dataKCNQ2,sex=="M")$Z~subset(TDT_dataKCNQ2,sex=="M")$Geno,alternative = "two.sided", paired=F)

Z.KCNQ40 <- lm(Z ~ Replica+Geno, data = subset(TDT_dataKCNQ40,sex=="F"))
summary(Z.KCNQ40)
anova(Z.KCNQ40)
wilcox.test (subset(TDT_dataKCNQ40,sex=="F")$Z~subset(TDT_dataKCNQ40,sex=="F")$Geno,alternative = "two.sided", paired=F)
wilcox.test (subset(TDT_dataKCNQ40,sex=="M")$Z~subset(TDT_dataKCNQ40,sex=="M")$Geno,alternative = "two.sided", paired=F)

# LMM ON KNOCKDOWN TIME PER SEX -----------------------------------
# performs all lmm and save them in a text file called LMER_KO_PER_SEX.txt

sink(file = "LMER_KO_PER_SEX.txt")
Genes <- list(datarobo3,datamam,datashot,dataKCNQ2)
nombresGenes <- c("robo3","mam","shot","KCNQ")
temperatures <- c("37","38","39","40")
for (g in 1:length(Genes)){
  for (s in c("F","M")){
    for (t in temperatures){
      cat("\n############################################\n")
      cat("GENE: ",nombresGenes[g],"   SEX: ",s, "   TEMPERATURE: ",t)
      cat("\n############################################\n")
      modelo <- lmer(KO ~ Geno+(1|Replica), data = subset(Genes[[g]],temp==t&sex==s))
      cat("\n######################## summary\n")
      print(summary(modelo))
      cat("\n######################## ANOVA\n")
      print(anova(modelo))
      cat("\n######################## RANDOM EFFECTS\n")
      print(rand(modelo))
      cat("\n######################## SHAPIRO TEST\n")
      print(shapiro.test(resid(modelo)))
      cat("\n######################## FLIGNER TEST\n")
      print(fligner.test(KO ~ Geno, data = subset(Genes[[g]],temp==t&sex==s)))
    }
  }
}
sink(file = NULL)


# BOX PLOTS--------------
# prepare data
datarobo3$Geno = str_replace_all(datarobo3$Geno,"CONTROL-attp40","Control")
datarobo3$Geno = str_replace_all(datarobo3$Geno,"robo3-attp40","RNAi")
datamam$Geno = str_replace_all(datamam$Geno,"CONTROL-kk","Control")
datamam$Geno = str_replace_all(datamam$Geno,"mam-kk","RNAi")
dataKCNQ2$Geno = str_replace_all(dataKCNQ2$Geno,"CONTROL-attp2","Control")
dataKCNQ2$Geno = str_replace_all(dataKCNQ2$Geno,"KCNQ-attp2","RNAi")
datashot$Geno = str_replace_all(dataKCNQ2$Geno,"CONTROL-sh","Control")
datashot$Geno = str_replace_all(dataKCNQ2$Geno,"shot-sh","RNAi")

# Prepare legend
forlegend <- ggboxplot(subset(datamam,sex=="F"), x="temp", y="KO",fill="Geno",palette="Dark2", xlab="Temperature (°C)", 
                       ylab="Knockdown time (min)", ylim=c(15,125),
                       legend="right")+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=15),
        legend.position = "bottom")

leyenda <- get_legend(forlegend)

## MAM
mamFemales <-ggboxplot(subset(datamam,sex=="F"), x="temp", y="KO", fill="Geno", palette="Dark2", xlab="Temperature (°C)", 
                       ylab="Knockdown time (min)", ylim=c(15,125),
                       legend="right")+
  labs(fill="mam")+
  theme(legend.title=element_text(size=20), legend.text = element_text(size=20),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(126, 126)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=129),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(67, 67)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=70),
            aes(x=x, y=y, label="*"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(35, 35)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=38),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(30, 30)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=33),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="mam"), size=5,
            inherit.aes=FALSE)

mamMales <-ggboxplot(subset(datamam,sex=="M"), x="temp", y="KO", fill="Geno", palette="Dark2", xlab="Temperature (°C)", 
                     ylab="Knockdown time (min)", ylim=c(15,125),
                     legend="right")+
  labs(fill="mam")+
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(120, 120)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=123),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(75, 75)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=78),
            aes(x=x, y=y, label="***"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(45, 45)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=48),
            aes(x=x, y=y, label="**"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(30, 30)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=33),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="mam"), size=5,
            inherit.aes=FALSE)

mam_final <- grid.arrange(mamFemales,mamMales,ncol=2, widths=c(1.2,1))


# ROBO3
robo3Female <- ggboxplot(subset(datarobo3,sex=="F"), x="temp", y="KO",fill="Geno",palette="Dark2", xlab="Temperature (°C)", 
                         ylab="Knockdown time (min)", ylim=c(15,125),
                         legend="right")+
  labs(fill="robo3")+
  theme(legend.title=element_text(size=20),legend.text = element_text(size=15),legend.position = "none")+ # erase legend from females
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(123, 123)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=126),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(75, 75)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=78),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(45, 45)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=48),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(30, 30)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=33),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="robo3"), size=5,
            inherit.aes=FALSE)

robo3male <- ggboxplot(subset(datarobo3,sex=="M"), x="temp", y="KO",fill="Geno",palette="Dark2", xlab="Temperature (°C)", 
                       ylab="Knockdown time (min)", ylim=c(15,125),
                       legend="right")+
  labs(fill="robo3")+
  theme(legend.title=element_text(size=20),legend.text = element_text(size=15),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(120, 120)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=123),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(75, 75)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=78),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(45, 45)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=48),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(30, 30)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=33),
            aes(x=x, y=y, label="*"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="robo3"), size=5,
            inherit.aes=FALSE)

robo3_final <- grid.arrange(robo3Female,robo3male, ncol=2, widths=c(1.2,1))

## KCNQ
KCNQFemales <-ggboxplot(subset(dataKCNQ2,sex=="F"), x="temp", y="KO", fill="Geno", palette="Dark2", xlab="Temperature (°C)", 
                        ylab="Knockdown time (min)", ylim=c(15,125),
                        legend="right")+
  labs(fill="KCNQ")+
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(127, 127)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=130),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(67, 67)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=70),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(40, 40)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=43),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(28, 28)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=31),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="KCNQ"), size=5,
            inherit.aes=FALSE)

KCNQMales <-ggboxplot(subset(dataKCNQ2,sex=="M"), x="temp", y="KO", fill="Geno", palette="Dark2", xlab="Temperature (°C)", 
                      ylab="Knockdown time (min)", ylim=c(15,125),
                      legend="right")+
  labs(fill="KCNQ")+
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(110, 110)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=113),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(77, 77)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=80),
            aes(x=x, y=y, label="***"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(47, 47)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=50),
            aes(x=x, y=y, label="**"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=34),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="KCNQ"), size=5,
            inherit.aes=FALSE)

KCNQ_final <- grid.arrange(KCNQFemales,KCNQMales, ncol=2, widths=c(1.2,1))

## SHOT
shotFemales <-ggboxplot(subset(datashot,sex=="F"), x="temp", y="KO", fill="Geno", palette="Dark2", xlab="Temperature (°C)", 
                        ylab="Knockdown time (min)", ylim=c(15,130),
                        legend="right")+
  labs(fill="shot")+
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank())+
  scale_y_continuous(breaks = seq(25, 125, 25))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(130, 130)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=133),
            aes(x=x, y=y, label="*"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(77, 77)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=80),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(47, 47)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=50),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=34),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=130), 
            aes(x=x, y=y, label="shot"), size=5,
            inherit.aes=FALSE)

shotMales <-ggboxplot(subset(datashot,sex=="M"), x="temp", y="KO", fill="Geno", palette="Dark2", xlab="Temperature (°C)", 
                      ylab="Knockdown time (min)", ylim=c(15,130),
                      legend="right")+
  labs(fill="shot")+
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank())+
  scale_y_continuous(breaks = seq(25, 125, 25))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(125, 125)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=128),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(76, 76)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=79),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(46, 46)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=49),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=34),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=130), 
            aes(x=x, y=y, label="shot"), size=5,
            inherit.aes=FALSE)

shot_final <- grid.arrange(shotFemales,shotMales, ncol=2, widths=c(1.2,1))

# Combine graphs
lay <- rbind(c(1,2),
             c(3,3),
             c(4,4),
             c(5,5),
             c(6,6),
             c(7,7))

# combine graphs and save them in png file called rnai_all.png
png(file="rnai_all.png",width=2200, height=3000, res = 300)
grid.arrange(grobs=list(text_grob("Females",size=15),text_grob("Males",size=15),mam_final, KCNQ_final,robo3_final,shot_final,leyenda), 
             layout_matrix=lay,left=text_grob("Knockdown time (min)", size=15,rot=90), heights=c(1,4,4,4,4.5,1))
dev.off()
