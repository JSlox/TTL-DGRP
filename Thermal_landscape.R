###################################################
# CALCULATE THERMAL TOLERANCE LANDSCAPE PARAMETERS
###################################################

# Function to calculate CTmax and z:

# TDT(matriX) -> data.frame with means of KO by temperature, sex, lines, CTmax y z
# receive a data.frame with column "KO" with KO times, "sex" with sex, "temp" with temperature, "DGRP" with DGRP line.
TDT <- function(matriz){
  #### CTmax and Z estimation ####
  
  # Calculate the knockdown time mean for each DGRP line and sex
  rep.st37 <- round(with(subset(matriz,temp=="37"),tapply(KO,list(DGRP,sex),mean)),2)
  rep.st38 <- round(with(subset(matriz,temp=="38"),tapply(KO,list(DGRP,sex),mean)),2)
  rep.st39 <- round(with(subset(matriz,temp=="39"),tapply(KO,list(DGRP,sex),mean)),2)
  rep.st40 <- round(with(subset(matriz,temp=="40"),tapply(KO,list(DGRP,sex),mean)),2)
  
  mean.rep <- with(matriz,tapply(KO,list(temp,DGRP,sex),mean))		# Calculate the knockdown time mean for each temperature, DGRP line and sex
  nlines <- length(unique(matriz$DGRP)) # extract number of DGRP lines from data set
  out <- matrix(,(nlines*2),6)	# Create an empty matrix with nlinesX2 rows (nlines x 2 sexes) and 6 columns (4 temps, CTmax and Z)
  for(i in 1:(nlines*2)){out[i,1] <- rep.st37[(i)]}	# Loop to paste the knockdown time means at 37ºC in the first column
  for(i in 1:(nlines*2)){out[i,2] <- rep.st38[(i)]}	# Loop to paste the knockdown time means at 38ºC in the second column
  for(i in 1:(nlines*2)){out[i,3] <- rep.st39[(i)]}	# Loop to paste the knockdown time means at 39ºC in the third column
  for(i in 1:(nlines*2)){out[i,4] <- rep.st40[(i)]}	# Loop to paste the knockdown time means at 40ºC in the fourth column
  for(i in 1:(nlines*2)){
    fit <- coef(lm(log10(mean.rep[(i*4-3):(i*4)]) ~ c(37:40)))	# Loop to calculate the intercept ans slope of each line
    fit <- c(-fit[1]/fit[2],1/fit[2])	# Transform intercept and slopr into CTmax and Z, respectively
    out[i,5] <- round(fit[1],2); out[i,6] <- round(-fit[2],2)} # Loop to move CTmax and Z values to fifth and sixth columns
  
  colnames(out) <- c("T37","T38","T39","T40","CTmax","Z")	# Assign column names
  sex <- factor(rep(c("F","M"),each=nlines))		# Create a factor for females and males
  DGRP <- c(levels(matriz$DGRP))	# Take DGRP levels and..
  DGRP <- factor(rep(DGRP,2))	# Create a factor with them
  tdt <- data.frame(DGRP,sex,out)
  return(tdt)
}