###################################
# THERMAL TOLERANCE LANDSCAPE PLOTS
###################################

#LOAD REQUIRED PACKAGES
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(readxl)

# READ DATA -------------------------------------------------------------------
data <- read_excel("KO_4temp_100DGRP.xlsx")

data$Geno <- as.factor(data$Geno)
data$temp <- as.factor(data$temp)
data$sex <- as.factor(data$sex)
data$Bloque <- as.factor(data$Bloque)
data$Fecha <- as.factor(data$Fecha)
data$Celda <- as.factor(data$Celda)
data$KO <- as.numeric(data$KO/60)

# CALCULATE CTMAX AND Z FROM ORIGINAL DATA ------------------------------------
source("Thermal_landscape.R")
TDT_table <- TDT(data)
TDT_table$sex <- as.factor(TDT_table$sex)

# MAKE TABLE WITH MEAN OF KO TIME BY LINE, SEX AND TEMPERATURE ----------------
rep.statics <- 10^with(data,tapply(log10(KO),list(temp,Geno,sex),mean))	#geometric mean for each replicate
out <- matrix(,800,1)
for(i in 1:800){out[i,1] <- rep.statics[(i)]}
colnames(out) <- c("TimeKO")
sex <- factor(rep(c(1,2),each=400),labels=c("females","males"))
Geno <- c(levels(data$Geno))
Geno <- factor(rep(c(1:100),each=4,2),labels=Geno)
temp <- factor(rep(c(37:40),200))
statics.rep <- data.frame(temp,Geno,sex,out)
data.F <- statics.rep[which(statics.rep[,"sex"]=="females"),]
data.M <- statics.rep[which(statics.rep[,"sex"]=="males"),]

# TDT PLOT --------------------------------------------------------------------
# Females
par(mfrow=c(1,2))
par(mar = c(5, 5, 1, 0) + 0.1)
plot(c(1, 4.05), c(4, 250),xlab="",ylab="",axes = FALSE,log="y",las=1,ylim=c(4,250),xlim=c(36.5,40.5),xaxs="i",yaxs="i")
for(i in 1:100){
  abline(lm(log10(data.F[(i*4-3):(i*4),4]) ~ c(37:40)),col="firebrick1", lwd=1)} 

axis(side = 1, at = c(37,38,39,40), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(5,10,20,50,100,200), cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Knockdown time (min)", cex.lab = 2, mgp = c(3.5, 0, 0))
title(xlab = "Temperature (°C)", cex.lab = 2, mgp = c(3.5, 0, 0))
legend("topright", bty = "n",legend = "Females", text.col = "firebrick1", x.intersp = 0.5, y.intersp = 1, cex = 2)

# Males
par(mar = c(5, 2, 1, 3) + 0.1)
plot(c(1, 4.05), c(4, 250),xlab="",ylab="",axes = FALSE,log="y",las=1,ylim=c(4,250),xlim=c(36.5,40.5),xaxs="i",yaxs="i")
for(i in 1:100){
  abline(lm(log10(data.M[(i*4-3):(i*4),4]) ~ c(37:40)),col="dodgerblue1", lwd=1)} 

axis(side = 1, at = c(37,38,39,40), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(5,10,20,50,100,200),labels=F, cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(xlab = "Temperature (°C)", cex.lab = 2, mgp = c(3.5, 0, 0))
legend("topright", bty = "n",legend = "Males", text.col = "dodgerblue1", x.intersp = 0.5, y.intersp = 1, cex = 2)
