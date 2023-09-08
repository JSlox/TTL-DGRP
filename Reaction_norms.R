###########3###################################
# Make reaction norms by temperatures and sexes
###############################################

# Read data ----------------
setwd("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/Github_TTL_DGRP/Data")
data <- read_excel("KO_4temp_100DGRP.xlsx")

data$DGRP <- as.factor(data$DGRP)
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
DGRP37.female <- tapply(subset(data37, sex == "F")$KO, subset(data37, sex == "F")$DGRP,mean)
DGRP37.male <- tapply(subset(data37, sex == "M")$KO, subset(data37, sex == "M")$DGRP,mean)
DGRP38.female <- tapply(subset(data38, sex == "F")$KO, subset(data38, sex == "F")$DGRP,mean)
DGRP38.male <- tapply(subset(data38, sex == "M")$KO, subset(data38, sex == "M")$DGRP,mean)
DGRP39.female <- tapply(subset(data39, sex == "F")$KO, subset(data39, sex == "F")$DGRP,mean)
DGRP39.male <- tapply(subset(data39, sex == "M")$KO, subset(data39, sex == "M")$DGRP,mean)
DGRP40.female <- tapply(subset(data40, sex == "F")$KO, subset(data40, sex == "F")$DGRP,mean)
DGRP40.male <- tapply(subset(data40, sex == "M")$KO, subset(data40, sex == "M")$DGRP,mean)

# REACTION NORMS-----------------------

# Females
par(mfrow=c(1,2))
par(mar = c(5, 4, 1, 2) + 0.1)
plot(c(1, 4.05), c(10, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in 1:100) {
  points(c(1, 2, 3, 4), c(DGRP37.female[i], DGRP38.female[i], DGRP39.female[i], DGRP40.female[i]), type = "b", cex = 1, lwd=1 , pch = 21, col = "firebrick1")
  points(c(1, 2, 3, 4), c(DGRP37.female[i], DGRP38.female[i], DGRP39.female[i], DGRP40.female[i]), type = "p", cex = 1.2, pch = 21, bg = brewer.pal(9, "Set1")[c(3, 1, 2,9)], col = NA)
}

axis(side = 1, at = c(1, 2, 3,4), labels = c("37","38","39","40"), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(10,45,80,115,150), cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Knockdown time (min)", cex.lab = 2, mgp = c(2.5, 0, 0))
title(xlab = "Temperature (°C)", cex.lab = 2, mgp = c(3.5, 0, 0))
legend("topright", bty = "n",legend = "Females", text.col = "firebrick1", x.intersp = 0.5, y.intersp = 1, cex = 2)


# Males
par(mar = c(5, 1, 1, 2) + 0.1)
plot(c(1, 4.05), c(10, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in 1:100) {
  points(c(1, 2, 3, 4), c(DGRP37.male[i], DGRP38.male[i], DGRP39.male[i], DGRP40.male[i]), type = "b", cex = 1, lwd=1, pch = 21, col = "dodgerblue1")
  points(c(1, 2, 3, 4), c(DGRP37.male[i], DGRP38.male[i], DGRP39.male[i], DGRP40.male[i]), type = "p", cex = 1.2, pch = 21, bg = brewer.pal(9, "Set1")[c(3, 1, 2,9)], col = NA)
}

axis(side = 1, at = c(1, 2, 3,4), labels = c("37","38","39","40"), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(10,45,80,115,150), labels=F, cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(xlab = "Temperature (°C)", cex.lab = 2, mgp = c(3.5, 0, 0))
legend("topright", bty = "n",legend = "Males", text.col = "dodgerblue1", x.intersp = 0.5, y.intersp = 1, cex = 2)


# sex by temperature
par(mar = c(5, 4, 1, 2) + 0.1)
plot(c(0.7, 4.3), c(10, 150), type = "n", axes = FALSE, xlab = "", ylab = "")

for (i in 1:100) {
  points(c(0.7, 1.3), c(DGRP37.female[i], DGRP37.male[i]), type = "b", cex = 1, pch = 21, col = brewer.pal(9, "Set1")[3])
  points(c(0.7, 1.3), c(DGRP37.female[i], DGRP37.male[i]), type = "p", cex = 1.2, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = NA)
}

for (i in 1:100) {
  points(c(1.7, 2.3), c(DGRP38.female[i], DGRP38.male[i]), type = "b", cex = 1, pch = 21, col = brewer.pal(9, "Set1")[1])
  points(c(1.7, 2.3), c(DGRP38.female[i], DGRP38.male[i]), type = "p", cex = 1.2, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = NA)
}

for (i in 1:100) {
  points(c(2.7, 3.3), c(DGRP39.female[i], DGRP39.male[i]), type = "b", cex = 1, pch = 21, col = brewer.pal(9, "Set1")[2])
  points(c(2.7, 3.3), c(DGRP39.female[i], DGRP39.male[i]), type = "p", cex = 1.2, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = NA)
}

for (i in 1:100) {
  points(c(3.7, 4.3), c(DGRP40.female[i], DGRP40.male[i]), type = "b", cex = 1, pch = 21, col = brewer.pal(9, "Set1")[9])
  points(c(3.7, 4.3), c(DGRP40.female[i], DGRP40.male[i]), type = "p", cex = 1.2, pch = 21, bg = c("firebrick1", "dodgerblue1"), col = NA)
}

axis(side = 1, at = c(1, 2, 3,4), labels = c("37","38","39","40"), cex.axis = 2, mgp = c(3, 1, 0), lwd = 1)
axis(side = 2, at = c(10,45,80,115,150), cex.axis = 2, mgp = c(3, 0.5, 0), lwd = 1)
box(bty = "l")

title(ylab = "Knockdown time (min)", cex.lab = 2, mgp = c(2.5, 0, 0))
title(xlab = "Temperature (°C)", cex.lab = 2, mgp = c(3.5, 0, 0))
legend("topright", bty = "n", pch = 21, col = c("firebrick1", "dodgerblue1"), pt.cex = 1.5, pt.bg = c("firebrick1", "dodgerblue1"), legend = c("Females", "Males"), text.col = c("firebrick1", "dodgerblue1"), x.intersp = 0.5, y.intersp = 1, cex = 2.2)