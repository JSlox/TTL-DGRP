################
# RNAi ANALYSIS 
################

#LOAD REQUIRED PACKAGES
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(readxl, reshape2, ggplot2, stringr, ggpubr, tidyverse, gridExtra, lme4, lmerTest, RColorBrewer)

# READ DATA -------------------------------------------------------------------
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
source("Thermal_landscape.R")

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

# LMM ON CTMAX AND Z ------------------------------------------

#CTmax y Z

TDT_datarobo3 <- droplevels(subset(TDT_full,Geno=="robo3-attp40"|Geno=="CONTROL-attp40"))
TDT_datamam <- droplevels(subset(TDT_full,Geno=="mam-kk"|Geno=="CONTROL-kk"))
TDT_datashot <- droplevels(subset(TDT_full,Geno=="shot-sh"|Geno=="CONTROL-sh"))
TDT_dataKCNQ2 <- droplevels(subset(TDT_full,Geno=="KCNQ-attp2"|Geno=="CONTROL-attp2"))

# CTmax
# LMM ON CTMAX FULL MODEL -----------------------------------
# performs all lmm and save them in a text file called LMER_CTMAX_FULL.txt

sink(file = "LMER_CTMAX_FULL.txt")
Genes <- list(TDT_datarobo3,TDT_datamam,TDT_datashot,TDT_dataKCNQ2)
nombresGenes <- c("robo3","mam","shot","KCNQ")
for (g in 1:length(Genes)){
    cat("\n############################################\n")
    cat("CTMAX   GENE: ",nombresGenes[g])
    cat("\n############################################\n")
    modelo <- lmer(CTmax ~ Geno*sex+(1|Replica), data = Genes[[g]])
    cat("\n######################## summary\n")
    print(summary(modelo))
    cat("\n######################## ANOVA\n")
    print(anova(modelo))
    cat("\n######################## RANDOM EFFECTS\n")
    print(rand(modelo))
    cat("\n######################## SHAPIRO TEST\n")
    print(shapiro.test(resid(modelo)))
    cat("\n######################## FLIGNER TEST\n")
    print(fligner.test(CTmax ~ Geno, data = Genes[[g]]))
}
sink(file = NULL)

# LMM ON CTmax PER SEX -----------------------------------
# performs all lmm and save them in a text file called LMER_CTmax_PER_SEX.txt
sink(file = "LMER_CTmax_PER_SEX.txt")
Genes <- list(TDT_datarobo3,TDT_datamam,TDT_datashot,TDT_dataKCNQ2)
nombresGenes <- c("robo3","mam","shot","KCNQ")
for (g in 1:length(Genes)){
  for (s in c("F","M")){
    cat("\n############################################\n")
    cat("GENE: ",nombresGenes[g],"   SEX: ",s)
    cat("\n############################################\n")
    modelo <- lmer(CTmax ~ Geno+(1|Replica), data = subset(Genes[[g]],sex==s))
    cat("\n######################## summary\n")
    print(summary(modelo))
    cat("\n######################## ANOVA\n")
    print(anova(modelo))
    cat("\n######################## RANDOM EFFECTS\n")
    print(rand(modelo))
    cat("\n######################## SHAPIRO TEST\n")
    print(shapiro.test(resid(modelo)))
    cat("\n######################## FLIGNER TEST\n")
    print(fligner.test(CTmax ~ Geno, data = subset(Genes[[g]],sex==s)))
  }
}
sink(file = NULL)

# Z
# LMM ON Z FULL MODEL -----------------------------------
# performs all lmm and save them in a text file called LMER_Z_FULL.txt

sink(file = "LMER_z_FULL.txt")
Genes <- list(TDT_datarobo3,TDT_datamam,TDT_datashot,TDT_dataKCNQ2)
nombresGenes <- c("robo3","mam","shot","KCNQ")
for (g in 1:length(Genes)){
  cat("\n############################################\n")
  cat("Z   GENE: ",nombresGenes[g])
  cat("\n############################################\n")
  modelo <- lmer(Z ~ Geno*sex+(1|Replica), data = Genes[[g]])
  cat("\n######################## summary\n")
  print(summary(modelo))
  cat("\n######################## ANOVA\n")
  print(anova(modelo))
  cat("\n######################## RANDOM EFFECTS\n")
  print(rand(modelo))
  cat("\n######################## SHAPIRO TEST\n")
  print(shapiro.test(resid(modelo)))
  cat("\n######################## FLIGNER TEST\n")
  print(fligner.test(Z ~ Geno, data = Genes[[g]]))
}
sink(file = NULL)

# LMM ON z PER SEX -----------------------------------
# performs all lmm and save them in a text file called LMER_z_PER_SEX.txt
sink(file = "LMER_z_PER_SEX.txt")
Genes <- list(TDT_datarobo3,TDT_datamam,TDT_datashot,TDT_dataKCNQ2)
nombresGenes <- c("robo3","mam","shot","KCNQ")
for (g in 1:length(Genes)){
  for (s in c("F","M")){
    cat("\n############################################\n")
    cat("GENE: ",nombresGenes[g],"   SEX: ",s)
    cat("\n############################################\n")
    modelo <- lmer(Z ~ Geno+(1|Replica), data = subset(Genes[[g]],sex==s))
    cat("\n######################## summary\n")
    print(summary(modelo))
    cat("\n######################## ANOVA\n")
    print(anova(modelo))
    cat("\n######################## RANDOM EFFECTS\n")
    print(rand(modelo))
    cat("\n######################## SHAPIRO TEST\n")
    print(shapiro.test(resid(modelo)))
    cat("\n######################## FLIGNER TEST\n")
    print(fligner.test(Z ~ Geno, data = subset(Genes[[g]],sex==s)))
  }
}
sink(file = NULL)

# LMM ON KNOCKDOWN TIME PER TEMPERATURE FULL MODEL -----------------------------------
# performs all lmm and save them in a text file called LMER_KO_PER_SEX.txt

sink(file = "LMER_KO_FULL.txt")
Genes <- list(datarobo3,datamam,datashot,dataKCNQ2)
nombresGenes <- c("robo3","mam","shot","KCNQ")
temperatures <- c("37","38","39","40")
for (g in 1:length(Genes)){
    for (t in temperatures){
      cat("\n############################################\n")
      cat("GENE: ",nombresGenes[g], "   TEMPERATURE: ",t)
      cat("\n############################################\n")
      modelo <- lmer(KO ~ Geno*sex+(1|Replica), data = subset(Genes[[g]],temp==t))
      cat("\n######################## summary\n")
      print(summary(modelo))
      cat("\n######################## ANOVA\n")
      print(anova(modelo))
      cat("\n######################## RANDOM EFFECTS\n")
      print(rand(modelo))
      cat("\n######################## SHAPIRO TEST\n")
      print(shapiro.test(resid(modelo)))
      cat("\n######################## FLIGNER TEST\n")
      print(fligner.test(KO ~ Geno, data = subset(Genes[[g]],temp==t)))
  }
}
sink(file = NULL)


# LMM ON KNOCKDOWN TIME PER TEMPERATURE AND SEX -----------------------------------
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

# CTMAX
# Prepare legend
forlegend <- ggboxplot(subset(datamam,sex=="F"), x="temp", y="KO",fill="Geno",palette="Dark2", xlab="Temperature (°C)", 
                       ylab="Knockdown time (min)", ylim=c(15,125),
                       legend="right")+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=15),
        legend.position = "bottom")

leyenda <- get_legend(forlegend)

## MAM
mamFemales <- ggplot(subset(datamam,sex=="F"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20), legend.text = element_text(size=20),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x= element_blank(),
        axis.text.y= element_text(colour="black"), axis.ticks = element_line(color = "black"))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(128, 128)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=132),
            aes(x=x, y=y, label="*"), size=3,
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
  geom_text(data=tibble(x=3, y=39),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(30, 30)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=34),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="mam"), size=5,
            inherit.aes=FALSE)

mamMales <-ggplot(subset(datamam,sex=="M"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks = element_line(color = "black"))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(122, 122)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=126),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(76, 76)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=79),
            aes(x=x, y=y, label="***"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(46, 46)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=49),
            aes(x=x, y=y, label="***"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(32, 32)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=36),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="mam"), size=5,
            inherit.aes=FALSE)

mam_final <- grid.arrange(mamFemales,mamMales,ncol=2, widths=c(1.2,1))


# ROBO3
robo3Female <- ggplot(subset(datarobo3,sex=="F"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20),legend.text = element_text(size=15),legend.position = "none")+ # erase legend from females
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15),axis.title.y = element_blank(),
        axis.text.y= element_text(colour="black"), axis.ticks = element_line(color = "black"),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(125, 125)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=129),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(77, 77)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=81),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(47, 47)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=51),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(32, 32)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=36),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="robo3"), size=5,
            inherit.aes=FALSE)

robo3male <- ggplot(subset(datarobo3,sex=="M"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20),legend.text = element_text(size=15),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks = element_line(color = "black"))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(120, 120)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=124),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(75, 75)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=79),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(47, 47)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=51),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=34),
            aes(x=x, y=y, label="**"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="robo3"), size=5,
            inherit.aes=FALSE)

robo3_final <- grid.arrange(robo3Female,robo3male, ncol=2, widths=c(1.2,1))

## KCNQ
KCNQFemales <-ggplot(subset(dataKCNQ2,sex=="F"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15), axis.title.y = element_blank(),
        axis.text.y= element_text(colour="black"), axis.ticks = element_line(color = "black"),
        axis.title.x = element_blank(),axis.text.x = element_blank())+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(129, 129)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=133),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(67, 67)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=71),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(40, 40)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=44),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(28, 28)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=32),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="KCNQ"), size=5,
            inherit.aes=FALSE)

KCNQMales <-ggplot(subset(dataKCNQ2,sex=="M"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks = element_line(color = "black"))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(110, 110)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=114),
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
            aes(x=x, y=y, label="***"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=35),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=125), 
            aes(x=x, y=y, label="KCNQ"), size=5,
            inherit.aes=FALSE)

KCNQ_final <- grid.arrange(KCNQFemales,KCNQMales, ncol=2, widths=c(1.2,1))

## SHOT
shotFemales <-ggplot(subset(datashot,sex=="F"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  xlab("Temperature (°C)")+
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135),breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15),
        axis.title.y = element_blank(),axis.text.y= element_text(colour="black"), axis.ticks = element_line(color = "black"),
        axis.text.x= element_text(colour="black"))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(131, 131)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=134),
            aes(x=x, y=y, label="*"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(77, 77)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=81),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(47, 47)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=51),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=35),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=130), 
            aes(x=x, y=y, label="shot"), size=5,
            inherit.aes=FALSE)

shotMales <-ggplot(subset(datashot,sex=="M"), aes(x=temp, y=KO, fill=Geno))+
  geom_boxplot(lwd=0.5, outlier.shape = NA, width=0.9) +
  theme_bw() +
  xlab("Temperature (°C)")+
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(limits= c(15,135), breaks = seq(25,125, by = 25))+
  geom_point(shape = 16,position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9), alpha = 0.4, size = 1) +
  theme(legend.title=element_text(size=20), legend.text = element_text(size=15),legend.position = "none")+
  theme(panel.border = element_blank(), axis.line = element_line(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size=15),axis.text = element_text(size=15), axis.text.y= element_blank(),
        axis.title.y = element_blank(),axis.text.x= element_text(colour="black"), axis.ticks = element_line(color = "black"))+
  geom_line(data=tibble(x=c(0.7, 1.3), y=c(125, 125)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=1, y=129),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(1.7, 2.3), y=c(76, 76)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=2, y=80),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2.7, 3.3), y=c(46, 46)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=3, y=50),
            aes(x=x, y=y, label="ns"), size=3,
            inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(3.7, 4.3), y=c(31, 31)),
            aes(x=x, y=y),
            inherit.aes=FALSE)+
  geom_text(data=tibble(x=4, y=35),
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

# combine graphs and save them in svg file called rnai_all.png
plot_grid <- grid.arrange(grobs=list(text_grob("Females",size=15),text_grob("Males",size=15),mam_final, KCNQ_final,robo3_final,shot_final,leyenda), 
             layout_matrix=lay,left=text_grob("Knockdown time (min)", size=15,rot=90), heights=c(1,4,4,4,4.5,1))
ggsave(file="rnai_all.svg", plot=plot_grid, width=2200, height=3000, dpi = 300, units="px")