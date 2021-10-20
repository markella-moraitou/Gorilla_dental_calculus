#Start logging
sink(file = "log8_alpha.txt")

#load packages
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(vegan)
library(MASS)
library(ggpubr)

load(".RData")

#Community-level taxonomic analysis - Script 8
#Comparing alpha diversity between samples

alpha <- estimate_richness(spe_data_final, measures = c("Observed", "Shannon", "Chao1"))

for (i in 1:nrow(alpha)){ #extract subspecies, seq. centre and read count metadata from spe_data_filt
  alpha$Spec.subspecies[i] = as.character(sample_data(spe_data_final)$Spec.subspecies[which(sample_names(spe_data_final)==row.names(alpha)[i])])
  alpha$Seq.centre[i] = as.character(sample_data(spe_data_final)$Seq.centre[which(sample_names(spe_data_final)==row.names(alpha)[i])])
  alpha$readcount[i] =  as.character(sample_data(spe_data_final)$readcount.m.before.Kraken[which(sample_names(spe_data_final)==row.names(alpha)[i])])
}
#Save plots
pdf(file = "/proj/sllstore2017021/nobackup/MARKELLA/T3_community-level/alpha_plots.pdf")

#Have a look at the data
ggboxplot(alpha, y = "Observed", x = "Spec.subspecies")
ggboxplot(alpha, y = "Shannon", x = "Spec.subspecies")
ggboxplot(alpha, y = "Chao1", x = "Spec.subspecies")
dev.off()

#### ANOVA ####
##For Observed richness

alpha$Spec.subspecies <- factor(alpha$Spec.subspecies, levels=c("gorilla", "graueri", "beringei"))
alpha$Seq.centre <- factor(alpha$Seq.centre)
alpha$readcount <- as.numeric(alpha$readcount)

#A square root transformation is needed
alpha$Observed_sqrt <- sqrt(500-alpha$Observed)

alphanova <- aov(data=alpha, Observed_sqrt ~ readcount + Seq.centre + Spec.subspecies)

#Check residuals
hist(studres(alphanova))
shapiro.test(studres(alphanova))

#Check for heteroskedasticity
plot(fitted.values(alphanova), studres(alphanova)) #There is no heteroscedasticity

summary(alphanova)

#Compare different levels using Tukey's test
TukeyHSD(alphanova, "Spec.subspecies")

##For Chao1 estimator

#A square root transformation is needed
alpha$Chao1_sqrt <- sqrt(500-alpha$Chao1)

chanova <- aov(data=alpha, Chao1_sqrt ~ readcount + Seq.centre + Spec.subspecies)

#Check residuals
hist(studres(chanova))
shapiro.test(studres(chanova))

#Check for heteroskedasticity
plot(fitted.values(chanova), studres(chanova))

summary(chanova)

#Compare different levels using Tukey's test
TukeyHSD(chanova, "Spec.subspecies")

#For Shannon index
#A exp transformation needed
alpha$Shannon_exp <- exp(alpha$Shannon)
shannova <- aov(data=alpha, Shannon_exp ~ readcount + Seq.centre + Spec.subspecies)

#Check residuals
hist(studres(shannova))
shapiro.test(studres(shannova))

#Check for heteroskedasticity
plot(fitted.values(shannova), studres(shannova)) 

summary(shannova)

#Compare different levels using Tukey's test
TukeyHSD(shannova, "Spec.subspecies")

sessionInfo()
#Stop logging
sink(file = NULL)

save.image()
