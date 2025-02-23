###############################################
# MAKE REACTION NORMS BY TEMPERATURES AND SEXES
###############################################

#LOAD REQUIRED PACKAGES
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(readxl, RColorBrewer)

# READ AND PREPARE DATA OF KO TIME -------------------------------------------------------------------
data <- read_excel("KO_4temp_100DGRP.xlsx")

data$Geno <- as.factor(data$Geno)
data$temp <- as.factor(data$temp)
data$sex <- as.factor(data$sex)
data$Bloque <- as.factor(data$Bloque)
data$Fecha <- as.factor(data$Fecha)
data$Celda <- as.factor(data$Celda)
data$KO <- as.numeric(data$KO/60)

data37 <- subset(data,temp=="37")
data38 <- subset(data,temp=="38")
data39 <- subset(data,temp=="39")
data40 <- subset(data,temp=="40")
DGRP37.female <- tapply(subset(data37, sex == "F")$KO, subset(data37, sex == "F")$Geno,mean)
DGRP37.male <- tapply(subset(data37, sex == "M")$KO, subset(data37, sex == "M")$Geno,mean)
DGRP38.female <- tapply(subset(data38, sex == "F")$KO, subset(data38, sex == "F")$Geno,mean)
DGRP38.male <- tapply(subset(data38, sex == "M")$KO, subset(data38, sex == "M")$Geno,mean)
DGRP39.female <- tapply(subset(data39, sex == "F")$KO, subset(data39, sex == "F")$Geno,mean)
DGRP39.male <- tapply(subset(data39, sex == "M")$KO, subset(data39, sex == "M")$Geno,mean)
DGRP40.female <- tapply(subset(data40, sex == "F")$KO, subset(data40, sex == "F")$Geno,mean)
DGRP40.male <- tapply(subset(data40, sex == "M")$KO, subset(data40, sex == "M")$Geno,mean)

# REACTION NORMS---------------------------------------------------------------


# prepares graph area
#par(mfrow=c(1,2)) 
tiff(file="Avg_reaction_norm_full.tiff",width=2000, height=2100,res=300,units = "px")
layout(matrix(c(1,1,2,3), nrow=2,byrow=TRUE),heights = c(1.5, 1))
par(mar = c(3.2, 3, 1, 2) + 0.1) #c(abajo, izquierda, arriba, derecha)

# sex by temperature -----------------------------------------------------------

# prepares graph area
#par(mfrow=c(1,1))
par(mar = c(3.2, 3, 1, 2) + 0.1)
# creates empty plot
plot(c(0.7, 4.3), c(10, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

# Graph lines and points for the both sexes at 37°C
for (i in 1:100) {
  points(c(0.7, 1.3), c(DGRP37.female[i], DGRP37.male[i]), type = "b", lwd=0.4 ,cex = 0.8, pch = 21, col = brewer.pal(9, "Set1")[3])
  points(c(0.7, 1.3), c(DGRP37.female[i], DGRP37.male[i]), type = "p", cex = 0.8, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = "black")
}

# Graph lines and points for the both sexes at 38°C
for (i in 1:100) {
  points(c(1.7, 2.3), c(DGRP38.female[i], DGRP38.male[i]), type = "b", lwd=0.4 ,cex = 0.8, pch = 21, col = brewer.pal(9, "Set1")[1])
  points(c(1.7, 2.3), c(DGRP38.female[i], DGRP38.male[i]), type = "p", cex = 0.8, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = "black")
}

# Graph lines and points for the both sexes at 39°C
for (i in 1:100) {
  points(c(2.7, 3.3), c(DGRP39.female[i], DGRP39.male[i]), type = "b", lwd=0.4 ,cex = 0.8, pch = 21, col = brewer.pal(9, "Set1")[2])
  points(c(2.7, 3.3), c(DGRP39.female[i], DGRP39.male[i]), type = "p", cex = 0.8, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = "black")
}

# Graph lines and points for the both sexes at 40°C
for (i in 1:100) {
  points(c(3.7, 4.3), c(DGRP40.female[i], DGRP40.male[i]), type = "b", lwd=0.4 ,cex = 0.8, pch = 21, col = brewer.pal(9, "Set1")[9])
  points(c(3.7, 4.3), c(DGRP40.female[i], DGRP40.male[i]), type = "p", cex = 0.8, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = "black")
}

# axis and labels
axis(side = 1, at = c(1, 2, 3,4), labels = c("37","38","39","40"), cex.axis = 1, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(10,45,80,115,150), cex.axis = 1, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Knockdown time (min)", cex.lab = 1.2, mgp = c(1.5, 0, 0))
title(xlab = "Temperature (°C)", cex.lab = 1.2, mgp = c(2.2, 0, 0))
legend("topright", bty = "n", pch = 21, col = c("firebrick1", "dodgerblue1"), pt.cex = 1.5, pt.bg = c("firebrick1", "dodgerblue1"), legend = c("Females", "Males"), text.col = c("firebrick1", "dodgerblue1"), x.intersp = 0.5, y.intersp = 1, cex = 1.5)


# Females
# creates empty plot
plot(c(1, 4.05), c(10, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

# Graph lines and points for the 4 temperatures
for (i in 1:100) {
  points(c(1, 2, 3, 4), c(DGRP37.female[i], DGRP38.female[i], DGRP39.female[i], DGRP40.female[i]), type = "b", cex = 0.8, lwd=0.4 , pch = 21, col = "firebrick1")
  points(c(1, 2, 3, 4), c(DGRP37.female[i], DGRP38.female[i], DGRP39.female[i], DGRP40.female[i]), type = "p", cex = 0.8, pch = 21, col="black", bg = brewer.pal(9, "Set1")[c(3, 1, 2,9)])
}

# axis and labels
axis(side = 1, at = c(1, 2, 3,4), labels = c("37","38","39","40"), cex.axis = 1, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(10,45,80,115,150), cex.axis = 1, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Knockdown time (min)", cex.lab = 1.2, mgp = c(1.5, 0, 0))
title(xlab = "Temperature (°C)", cex.lab = 1.2, mgp = c(2.2, 0, 0))
legend("topright", bty = "n",legend = "Females", text.col = "firebrick1", x.intersp = 0.5, y.intersp = 1, cex = 1.5)

# Males
par(mar = c(3.2, 1, 1, 2) + 0.1)
plot(c(1, 4.05), c(10, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in 1:100) {
  points(c(1, 2, 3, 4), c(DGRP37.male[i], DGRP38.male[i], DGRP39.male[i], DGRP40.male[i]), type = "b", cex = 0.8, lwd=0.4, pch = 21, col = "dodgerblue1")
  points(c(1, 2, 3, 4), c(DGRP37.male[i], DGRP38.male[i], DGRP39.male[i], DGRP40.male[i]), type = "p", cex = 0.8, pch = 21, bg = brewer.pal(9, "Set1")[c(3, 1, 2,9)], col = "black")
}

axis(side = 1, at = c(1, 2, 3,4), labels = c("37","38","39","40"), cex.axis = 1, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(10,45,80,115,150), labels=F, cex.axis = 1, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(xlab = "Temperature (°C)", cex.lab = 1.2, mgp = c(2.2, 0, 0))
legend("topright", bty = "n",legend = "Males", text.col = "dodgerblue1", x.intersp = 0.5, y.intersp = 1, cex = 1.5)

dev.off()
