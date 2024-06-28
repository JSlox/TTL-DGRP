###########################
# PLOT GENETIC CORRELATION
###########################

#LOAD REQUIRED PACKAGES
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(ggpubr)

# READ DATA -------------------------------------------------------------------
data <- read.table("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/Genetic_Correlation.txt", header=T)

data$line <- as.factor(data$line)
data$sex <- as.factor(data$sex)

## Phenotypic correlations
data2 <- data[,c(5:8)]
cor(data2)
cov(data2)

# mean per line and temperature
line.st37 <- c(with(data,tapply(st37,list(line),mean)))
line.st38 <- c(with(data,tapply(st38,list(line),mean)))
line.st39 <- c(with(data,tapply(st39,list(line),mean)))
line.st40 <- c(with(data,tapply(st40,list(line),mean)))

# Standard error of the mean
se.st37 <- c(with(data,tapply(st37,list(line),sd))/sqrt(with(data,tapply(st37,list(line),length))))
se.st38 <- c(with(data,tapply(st38,list(line),sd))/sqrt(with(data,tapply(st38,list(line),length))))
se.st39 <- c(with(data,tapply(st39,list(line),sd))/sqrt(with(data,tapply(st39,list(line),length))))
se.st40 <- c(with(data,tapply(st40,list(line),sd))/sqrt(with(data,tapply(st40,list(line),length))))

# new data frame with mean and standard error
mean.line <- data.frame(line.st37,se.st37,line.st38,se.st38,line.st39,se.st39,line.st40,se.st40)
mean.line2 <- mean.line[,c(1,3,5,7)]

## Genetic correlations
cor(mean.line2)
cov(mean.line2)

## Plot A 37 vs 38
plotA <- ggplot(data = mean.line,aes(x = line.st37,y = line.st38)) + 
  geom_errorbar(aes(ymin = line.st38-se.st38,ymax = line.st38+se.st38), color="red") + 
  geom_errorbarh(aes(xmin = line.st37-se.st37,xmax = line.st37+se.st37), color="red") +
  geom_point(size=2) +  theme_classic() +
  labs(x="Knockdown time at 37ºC (min)", y="Knockdown time at 38ºC (min)") +
  theme(axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"))
plotA
## Plot B 37 vs 39
plotB <- ggplot(data = mean.line,aes(x = line.st37,y = line.st39)) + 
  geom_errorbar(aes(ymin = line.st39-se.st39,ymax = line.st39+se.st39), color="red") + 
  geom_errorbarh(aes(xmin = line.st37-se.st37,xmax = line.st37+se.st37), color="red") +
  geom_point(size=2) +  theme_classic() +
  labs(x="Knockdown time at 37ºC (min)", y="Knockdown time at 39ºC (min)") +
  theme(axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"))
plotB
## Plot C 37 vs 40
plotC <- ggplot(data = mean.line,aes(x = line.st37,y = line.st40)) + 
  geom_errorbar(aes(ymin = line.st40-se.st40,ymax = line.st40+se.st40), color="red") + 
  geom_errorbarh(aes(xmin = line.st37-se.st37,xmax = line.st37+se.st37), color="red") +
  geom_point(size=2) +  theme_classic() +
  labs(x="Knockdown time at 37ºC (min)", y="Knockdown time at 40ºC (min)") +
  theme(axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"))
plotC
## Plot D 38 vs 39
plotD <- ggplot(data = mean.line,aes(x = line.st38,y = line.st39)) + 
  geom_errorbar(aes(ymin = line.st39-se.st39,ymax = line.st39+se.st39), color="red") + 
  geom_errorbarh(aes(xmin = line.st38-se.st38,xmax = line.st38+se.st38), color="red") +
  geom_point(size=2) +  theme_classic() +
  labs(x="Knockdown time at 38ºC (min)", y="Knockdown time at 39ºC (min)") +
  theme(axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"))
plotD
## Plot E 38 vs 40
plotE <- ggplot(data = mean.line,aes(x = line.st38,y = line.st40)) + 
  geom_errorbar(aes(ymin = line.st40-se.st40,ymax = line.st40+se.st40), color="red") + 
  geom_errorbarh(aes(xmin = line.st38-se.st38,xmax = line.st38+se.st38), color="red") +
  geom_point(size=2) +  theme_classic() +
  labs(x="Knockdown time at 38ºC (min)", y="Knockdown time at 40ºC (min)") +
  theme(axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"))
plotE
## Plot F 39 vs 40
plotF <- ggplot(data = mean.line,aes(x = line.st39,y = line.st40)) + 
  geom_errorbar(aes(ymin = line.st40-se.st40,ymax = line.st40+se.st40), color="red") + 
  geom_errorbarh(aes(xmin = line.st39-se.st39,xmax = line.st39+se.st39), color="red") +
  geom_point(size=2) +  theme_classic() +
  labs(x="Knockdown time at 39ºC (min)", y="Knockdown time at 40ºC (min)") +
  theme(axis.text.x = element_text(size=10, colour = "black"),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.title.x = element_text(size=12, colour = "black"),
        axis.title.y = element_text(size=12, colour = "black"))
plotF

# Arrange all plots
ggarrange(plotA,plotB,plotC,plotD,plotE,plotF, ncol=2, nrow=3,
          labels="AUTO", label.x=0.09)