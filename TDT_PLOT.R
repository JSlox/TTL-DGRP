###################################
# THERMAL TOLERANCE LANDSCAPE PLOTS
###################################

# READ DATA ----------------
library(readxl)
data <- read_excel("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/KO_4temp_100DGRP.xlsx")

data$DGRP <- as.factor(data$DGRP)
data$temp <- as.factor(data$temp)
data$sex <- as.factor(data$sex)
data$Bloque <- as.factor(data$Bloque)
data$Fecha <- as.factor(data$Fecha)
data$Celda <- as.factor(data$Celda)
data$KO <- as.numeric(data$KO/60)

# CALCULATE CTMAX AND Z FROM ORIGINAL DATA -------
source("Thermal_landscape.R")
TDT_table <- TDT(data)

summary(TDT_table)

##sexes comparison of CTmax and z
TDT_table$sex <- as.factor(TDT_table$sex)

#CTmax
shapiro.test(TDT_table$CTmax)
fligner.test(TDT_table$CTmax,TDT_table$sex)
wilcox.test(subset(TDT_table, sex=="F")$CTmax, subset(TDT_table, sex=="M")$CTmax, alternative = "two.sided")

#Z
shapiro.test(TDT_table$Z)
fligner.test(TDT_table$Z,TDT_table$sex)
wilcox.test(subset(TDT_table, sex=="F")$Z, subset(TDT_table, sex=="M")$Z, alternative = "two.sided")

cor(TDT_table$CTmax,TDT_table$Z)
anova(manova(cbind(CTmax,Z)~sex, data= TDT_table))#no se cumplen supuestos
summary(manova(cbind(CTmax,Z)~sex, data= TDT_table))# no se cumplen supuestos

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
rep.statics <- 10^with(data,tapply(log10(KO),list(temp,DGRP,sex),mean))	#geometric mean for each replicate
out <- matrix(,800,1)
for(i in 1:800){out[i,1] <- rep.statics[(i)]}
colnames(out) <- c("TimeKO")
sex <- factor(rep(c(1,2),each=400),labels=c("females","males"))
DGRP <- c(levels(data$DGRP))
DGRP <- factor(rep(c(1:100),each=4,2),labels=DGRP)
temp <- factor(rep(c(37:40),200))
statics.rep <- data.frame(temp,DGRP,sex,out)


### TDT curve plot
rep.statics <- 10^with(data,tapply(log10(KO),list(temp,DGRP,sex),mean))	#geometric mean for each replicate
data.F <- statics.rep[which(statics.rep[,"sex"]=="females"),]
data.M <- statics.rep[which(statics.rep[,"sex"]=="males"),]
head(data.F)
head(data.M)

# Females
png(file="TDTFyM.png",width=3000, height=2000, res = 300)
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
#dev.off()

# Males
#png(file="TDTmachos.png",width=3000, height=1700, res = 300)
par(mar = c(5, 2, 1, 3) + 0.1)
plot(c(1, 4.05), c(4, 250),xlab="",ylab="",axes = FALSE,log="y",las=1,ylim=c(4,250),xlim=c(36.5,40.5),xaxs="i",yaxs="i")
for(i in 1:100){
  abline(lm(log10(data.M[(i*4-3):(i*4),4]) ~ c(37:40)),col="dodgerblue1", lwd=1)} 

axis(side = 1, at = c(37,38,39,40), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(5,10,20,50,100,200),labels=F, cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

#title(ylab = "Knockdown time (min)", cex.lab = 2, mgp = c(2.5, 0, 0))
title(xlab = "Temperature (°C)", cex.lab = 2, mgp = c(3.5, 0, 0))
legend("topright", bty = "n",legend = "Males", text.col = "dodgerblue1", x.intersp = 0.5, y.intersp = 1, cex = 2)
dev.off()