######################
# LINEAL MIXED MODELS
######################

#LOAD REQUIRED PACKAGES
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(readxl, lme4, lmerTest)

# READ DATA -------------------------------------------------------------------
data <- read_excel("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/scripts/Data/KO_4temp_100DGRP.xlsx")

data$Geno <- as.factor(data$Geno)
data$temp <- as.factor(data$temp)
data$sex <- as.factor(data$sex)
data$Bloque <- as.factor(data$Bloque)
data$Fecha <- as.factor(data$Fecha)
data$Celda <- as.factor(data$Celda)
data$KO <- as.numeric(data$KO/60)

# LMMs ------------------------------------------------------------------------

full.model <- lmer(KO ~ temp + sex + Bloque +temp:sex + (1|Geno) +(1|Geno:temp) + (1|Geno:sex) + (1|Geno:temp:sex), data = data)
summary(full.model)
anova(full.model)
rand(full.model)

# Model per sex
females.model <- lmer(KO ~ temp +Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex=="F"))
summary(females.model, correlation=F)
rand(females.model)
anova(females.model)

males.model <- lmer(KO ~ temp +Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex=="M"))
summary(males.model, correlation=F)
rand(males.model)
anova(males.model)

# 37°C models
temp37.model <- lmer(KO ~ sex +Bloque + (1|Geno) + (1|Geno:sex), data = subset(data,temp=="37"))
summary(temp37.model)
rand(temp37.model)
anova(temp37.model)

temp37females.model <- lmer(KO ~ Bloque +(1|Geno), data = subset(subset(data,temp=="37"),sex=="F"))
summary(temp37females.model)
rand(temp37females.model)
anova(temp37females.model)

temp37males.model <- lmer(KO ~ Bloque +(1|Geno), data = subset(subset(data,temp=="37"),sex=="M"))
summary(temp37males.model)
rand(temp37males.model)
anova(temp37males.model)

# 38°C models
temp38.model <- lmer(KO ~ sex +Bloque +(1|Geno) + (1|Geno:sex), data = subset(data,temp=="38"))
summary(temp38.model)
rand(temp38.model)
anova(temp38.model)

temp38females.model <- lmer(KO ~ Bloque+(1|Geno), data = subset(subset(data,temp=="38"),sex=="F"))
summary(temp38females.model)
rand(temp38females.model)
anova(temp38females.model)

temp38males.model <- lmer(KO ~ Bloque+(1|Geno), data = subset(subset(data,temp=="38"),sex=="M"))
summary(temp38males.model)
rand(temp38males.model)
anova(temp38males.model)

# 39°C models
temp39.model <- lmer(KO ~ sex +Bloque+(1|Geno) + (1|Geno:sex), data = subset(data,temp=="39"))
summary(temp39.model)
rand(temp39.model)
anova(temp39.model)

temp39females.model <- lmer(KO ~ Bloque+(1|Geno), data = subset(subset(data,temp=="39"),sex=="F"))
summary(temp39females.model)
rand(temp39females.model)
anova(temp39females.model)

temp39males.model <- lmer(KO ~ Bloque+(1|Geno), data = subset(subset(data,temp=="39"),sex=="M"))
summary(temp39males.model)
rand(temp39males.model)
anova(temp39males.model)

# 40°C models
temp40.model <- lmer(KO ~ sex+Bloque+(1|Geno) + (1|Geno:sex), data = subset(data,temp=="40"))
summary(temp40.model)
rand(temp40.model)
anova(temp40.model)

temp40females.model <- lmer(KO ~ Bloque+(1|Geno), data = subset(subset(data,temp=="40"),sex=="F"))
summary(temp40females.model)
rand(temp40females.model)
anova(temp40females.model)

temp40males.model <- lmer(KO ~ Bloque+(1|Geno), data = subset(subset(data,temp=="40"),sex=="M"))
summary(temp40males.model)
rand(temp40males.model)
anova(temp40males.model)

# models for pairs of temperatures

# 37-38°C
female3738.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="F" & (temp=="37"| temp=="38")))
summary(female3738.model)
rand(female3738.model)
anova(female3738.model)

male3738.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="M" & (temp=="37"| temp=="38")))
summary(male3738.model)
rand(male3738.model)
anova(male3738.model)

# 37-39°C
female3739.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="F" & (temp=="37"| temp=="39")))
summary(female3739.model)
rand(female3739.model)
anova(female3739.model)

male3739.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="M" & (temp=="37"| temp=="39")))
summary(male3739.model)
rand(male3739.model)
anova(male3739.model)

# 37-40°C
female3740.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="F" & (temp=="37"| temp=="40")))
summary(female3740.model)
rand(female3740.model)
anova(female3740.model)

male3740.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="M" & (temp=="37"| temp=="40")))
summary(male3740.model)
rand(male3740.model)
anova(male3740.model)

# 38-39°C
female3839.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="F" & (temp=="38"| temp=="39")))
summary(female3839.model)
rand(female3839.model)
anova(female3839.model)

male3839.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="M" & (temp=="38"| temp=="39")))
summary(male3839.model)
rand(male3839.model)
anova(male3839.model)

# 38-40°C
female3840.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="F" & (temp=="38"| temp=="40")))
summary(female3840.model)
rand(female3840.model)
anova(female3840.model)

male3840.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="M" & (temp=="38"| temp=="40")))
summary(male3840.model)
rand(male3840.model)
anova(male3840.model)

# 39-40°C
female3940.model <- lmer(KO ~ temp +Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="F" & (temp=="39"| temp=="40")))
summary(female3940.model)
rand(female3940.model)
anova(female3940.model)

male3940.model <- lmer(KO ~ temp+Bloque+ (1|Geno) + (1|Geno:temp), data = subset(data,sex =="M" & (temp=="39"| temp=="40")))
summary(male3940.model)
rand(male3940.model)
anova(male3940.model)


# CALCULATE CTMAX AND Z FROM ORIGINAL DATA ------------------------------------
source("C:/Users/Juan/Documents/Postgrado/LAB_Bio_integrativa/Papers/Soto_etal_2022_TDT_DGRP/TTL-DGRP/Thermal_landscape.R")
TDT_table <- TDT(data)
TDT_table$sex <- as.factor(TDT_table$sex)

# CTmax
CTmax.model <- lmer(CTmax ~ sex + (1|Geno), data = TDT_table)
summary(CTmax.model)
rand(CTmax.model)
anova(CTmax.model)

# Z
z.model <- lmer(Z ~ sex + (1|Geno), data = TDT_table)
summary(z.model)
rand(z.model)
anova(z.model)