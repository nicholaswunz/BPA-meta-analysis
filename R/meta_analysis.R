# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 23/04/2020
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: BPA meta-analysis
# Description: Meta-analysis, meta-regressions and figure production

# Load packages
library(metafor) # for meta-analysis and meta-regressions
library(dplyr) # for data manipulation like pipes %>%
library(ggplot2)
library(cowplot)
mytheme <- theme_bw() + {theme(panel.border = element_blank(), # Border around plotting area.
                               panel.grid.major = element_blank(), # Major grid lines blank
                               panel.grid.minor = element_blank(), # Minor grid lines blank
                               axis.line = element_line(colour = "black", size = 0.8), # axis lin size
                               axis.ticks = element_line(colour = "black", size = 0.8),
                               axis.text = element_text(size = 10, colour = "black"), # axis text size
                               axis.title=element_text(size = 10), #axis title size
                               panel.background = element_rect(fill = "transparent"), # bg of the panel
                               plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                               legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                               legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg)
} # set up plot theme

# set directory
setwd('D:/Projects/Plastic pollution project - Sydney/Meta-analysis')

# load data
meta <- read.csv("meta data raw.csv")
str(meta) # check data strings

## Summary of data ##-------------------------------------------------------------------------
meta <- subset(meta, n_control != 1) # remove studies with sample size of 1
meta$lnRR_variance2 <- abs(meta$lnRR_variance) # change negative values to absolute values (for variance)
meta$lnRR_variance2[meta$lnRR_variance2 == 0] <- meta$lnRR_variance2 + 0.0001
meta <- subset(meta, first_author != "Alexander") # remove Alexander et al 1988 (showed time-lag bias)

str(meta)

nrow(meta) # number of effect sizes
length(unique(meta$study_ID)) # number of studies
length(unique(meta$scientific_name_OTL))# number of species

# how many species per taxa
meta %>% 
  group_by(taxa) %>% 
  summarise(count = n_distinct(scientific_name_OTL))

# ES counts for traits
meta %>% 
  group_by(traits) %>% 
  summarise(count = n_distinct(effect_size_ID))

# n of study ID by trait
meta %>% 
  group_by(trait) %>% 
  summarise(count = n_distinct(study_ID))

## Meta-analysis ##-------------------------------------------------------------------------
# Overall effect of BPA exposure (non-phylogenetic component)
overall.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2,
                       random = list(~1 | scientific_name, ~1 | study_ID, ~1 | effect_size_ID), 
                       method = "REML", data = meta)
summary(overall.m)

# calculate R2 when accounting for all moderators
fix <- var(as.numeric(as.vector(overall.m$b) %*% t(as.matrix(overall.m$X))))
R2m <- fix / (fix + sum(overall.m$sigma2))
100 * R2m # account 23.14% of varience

# Overall effect of BPA exposure (phylogenetic component)
# Construct phylogenetic tree from meta_phylo.R

# Loading phylogenetic matrix "phylo_cor"
load("Phylogeny/phylo_cor.Rdata") # phylo_cor

phylo.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2,
                  random = list(~1 | scientific_name, ~1 | scientific_name_OTL, ~1 | study_ID, ~1 | effect_size_ID), 
                  R = list(scientific_name_OTL = phylo_cor), method = "REML", data = meta)
summary(phylo.m)

# phylogenetic signal
library(phytools)
trait <- meta[,45]
names(trait) <- meta[,11]

# test Pagels lambda with 999 randomizations:
phylosig(phylo_branch, trait, method = "lambda", test = TRUE, nsim = 999)

AIC(overall.m, phylo.m)
# Phylogenetic signal is low, and AIC is lower in non-phylogenetic model
# Proceed with non-phylogenetic meta-regression

# obtain estimates
overall.m2 <- data.frame(estimate = overall.m$b, ci.lb = overall.m$ci.lb, ci.ub = overall.m$ci.ub)

## Overall I2 analysis ##-------------------------------------------------------------------------
names(overall.m) # check names of model output

# s^2_m: variation within study = random error
s2m <- sum(1/overall.m2$vi) * (overall.m2$k-1) / (sum(1/overall.m2$vi)^2 - sum((1/overall.m2$vi)^2)) 
# s^2_t = total variance
s2t <- sum(overall.m2$sigma2) + s2m 

sum(overall.m2$sigma2) / s2t # Total heterogeneity
overall.m2$sigma2[1] / s2t # Species = non-phylogenetic variance 
overall.m2$sigma2[2] / s2t # Group_ID = group (exp) variance 
overall.m2$sigma2[3] / s2t # ES_ID = effect size level variance 

## Publication bias and sensitivity analysis ##-------------------------------------------
# 1. Time-lag effect
library(ggpmisc)
ggplot(sublethal, aes(x = year_published, y = lnRR_direct, size = n_control)) + 
  geom_point(shape = 21, fill = "#4292c6", alpha = 0.5) + 
  labs(x = "Publication year", y = "Effect size", size = "N") +
  scale_size_area(max_size = 10) + mytheme +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  geom_hline(yintercept = 0, lty = 2) + 
  geom_smooth(method = "lm", size = 1, se = F, colour = "#084594") +
  stat_poly_eq(formula = y ~ x, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE)

pub.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ year_published, 
                       random = list(~1 | scientific_name, ~1 | study_ID, ~1 | effect_size_ID), 
                       method = "REML", data = meta)
summary(pub.m)

# 2. Funnel & trim-fill plot
resid <- rstandard(overall.m) # recover residuals, se, z, slab
resid.df <- do.call("rbind", resid)
resid.df2 <- t(resid.df)

# run an rma.uni on the residuals data and run either trimfill and regtest
resid.es <- rma(yi = resid, vi = se, data = resid.df2)
par(mfrow = c(1,2))
funnel(resid.es, yaxis = "seinv", legend = TRUE, main = "Inverse standard error")
taf <- trimfill(resid.es, estimator = "R0")
taf
funnel(taf, yaxis="seinv", legend=TRUE, main = "Trim-fill output")

# 3. Egger's regression
regtest(resid.es, model="lm")

# 4. sensitivity analysis
overall.cooks <- cooks.distance(overall.m)
dev.off()
plot(overall.cooks, type = "o", pch = 19, xlab = "Observed Outcome", ylab = "Cook's Distance")

## Trait overall ##------------------------------------------------------------------
# within traits alone
trait.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ trait-1, 
                   random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                   data = meta)
summary(trait.m)

# Maginal R2
fix <- var(as.numeric(as.vector(trait.m$b) %*% t(as.matrix(trait.m$X))))
R2m <- fix / (fix + sum(trait.m$sigma2))
100 * R2m

# output model and plot
trait.m2 <- data.frame(trait = substr(row.names(trait.m$b), 6, 100), #remove the word response
                        estimate = trait.m$b, ci.lb = trait.m$ci.lb, ci.ub = trait.m$ci.ub)
library(colortools)
wheel("#084594", num = 12) # generate 12 colour wheels from blue
traitcol <- wheel("#084594", num = 12)

ggplot(trait.m2, aes(x = trait, y = estimate, colour = trait)) +
  mytheme + geom_hline(yintercept=0, linetype = "dashed") + 
  geom_point(size = 4, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  xlab(NULL) + ylab("effect size (lnRR)") + theme(legend.position = "top") +
  scale_color_manual(values = traitcol) +
  theme(axis.text.x = element_text(size = 8, angle = 15, hjust = 1))

# Trait I2 analysis
names(overall.m) # check names of model output

# s^2_m: variation within study = random error
s2m <- sum(1/trait.m$vi) * (trait.m$k-1) / (sum(1/trait.m$vi)^2 - sum((1/trait.m$vi)^2)) 
# s^2_t = total variance
s2t <- sum(trait.m$sigma2) + s2m 

sum(trait.m$sigma2) / s2t # Total heterogeneity
trait.m$sigma2[1] / s2t # Species = non-phylogenetic varianc
trait.m$sigma2[2] / s2t # Group_ID = group (exp) variance
trait.m$sigma2[3] / s2t # ES_ID = effect size level variance

## Trait overall (excluding fish) ##------------------------------------------------------------------
# within traits alone
meta.nofish <- meta %>% filter(taxa != "Fish")

levels(meta.nofish$taxa)

meta.nofish %>% 
  group_by(trait) %>% 
  summarise(count = n_distinct(effect_size_ID))
             
nofish.trait.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ trait-1, 
                  random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                  data = meta.nofish)
summary(nofish.trait.m)

# output model and plot
nofish.trait.m2 <- data.frame(trait = substr(row.names(nofish.trait.m$b), 6, 100), #remove the word response
                       estimate = nofish.trait.m$b, ci.lb = nofish.trait.m$ci.lb, ci.ub = nofish.trait.m$ci.ub)

ggplot(nofish.trait.m2, aes(x = trait, y = estimate, colour = trait)) +
  mytheme + geom_hline(yintercept=0, linetype = "dashed") + 
  geom_point(size = 4, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  xlab(NULL) + ylab("effect size (lnRR)") + theme(legend.position = "top") +
  scale_color_manual(values = traitcol) +
  theme(axis.text.x = element_text(size = 8, angle = 15, hjust = 1))

## Trait sub-analyses ##--------------------------------------------------------
# 1. Abnormality ##
abnorm <- subset(meta, trait %in% "Abnormality & damage")

abnorm %>% # ES counts for abnormal responses
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

abnorm %>% # n of study ID by abnormal responses
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

# select 5 most well-sampled responses based on study ID
abnorm2 <- subset(abnorm, response %in% c("Deformity","Pericardial edema","ROS production", "DNA damage", "Swim bladder area"))

abnorm.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                  random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                  data = abnorm2)
summary(abnorm.m)

# extract model output and plot
abnorm.m2 <- data.frame(response = substr(row.names(abnorm.m$b), 9, 100), #remove the word response
                       estimate = abnorm.m$b, ci.lb = abnorm.m$ci.lb, ci.ub = abnorm.m$ci.ub)

abnorm.plot <- ggplot(abnorm.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept=0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE)+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 0.8, width=0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#084594","#396aa9", "#6b8fbf", "#9cb5d4", "#cedaea")) + #blue theme
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size=8, angle = 10, hjust = 1))

# 2. Behaviour ##
behav<- subset(meta, trait == "Behaviour")

behav %>% # ES counts for behaviour responses
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

behav %>% # n of study ID by behaviour responses
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

# select 5 most well-sampled responses based on study ID
behav2 <- subset(behav, response %in% c("Activity","Feeding","Distance moved", "Conspecific interaction", "Mate choice"))

behav.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                   random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                   data = behav2)
summary(behav.m)

# extract model output and plot
behav.m2 <- data.frame(response = substr(row.names(behav.m$b), 9, 100), #remove the word response
                        estimate = behav.m$b, ci.lb = behav.m$ci.lb, ci.ub = behav.m$ci.ub)

behav.plot <- ggplot(behav.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#110894", "#4139a9", "#706bbf", "#a09cd4", "#cfceea")) + # dark blue theme
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 3. Cardiovascular ##
cardio<- subset(meta, trait == "Cardiovascular")

cardio %>% # ES counts for cardio responses
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

cardio %>% # n of study ID by cardio responses
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

# select 5 most well-sampled responses based on study ID
cardio2 <- subset(cardio, response %in% c("Heart rate","Blood cell count","Haematocrit", "Haemoglobin", "Blood flow"))

cardio.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                  random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                  data = cardio2)
summary(cardio.m)

# extract model output and plot
cardio.m2 <- data.frame(response = substr(row.names(cardio.m$b), 9, 100), #remove the word response
                       estimate = cardio.m$b, ci.lb = cardio.m$ci.lb, ci.ub = cardio.m$ci.ub)

cardio.plot <- ggplot(cardio.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#570894", "#7939a9", "#9a6bbf", "#bc9cd4", "#ddceea")) + # purple theme
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 4. Development ##
dev <- subset(meta, trait == "Development")

dev %>% # ES counts for development reponses
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

dev %>% # n of study ID by development reponses
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

# select 5 most well-sampled responses based on study ID
dev2 <- subset(dev, response %in% c("Hatching success","Hatching time","Metamorphosis", "Moulting frequency", "Emergence"))

dev.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                   random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                   data = dev2)
summary(dev.m)

# extract model output and plot
dev.m2 <- data.frame(response = substr(row.names(dev.m$b), 9, 100), #remove the word response
                        estimate = dev.m$b, ci.lb = dev.m$ci.lb, ci.ub = dev.m$ci.ub)

dev.plot <- ggplot(dev.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#94088B","#a939a2", "#bf6bb9", "#d49cd1", "#eacee8")) + # magenta theme
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 5. Energy metabolism##
metab <- subset(meta, trait == "Energy metabolism")

metab %>% # ES counts for development reponses
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

metab %>% # n of study ID by development reponses
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

metab.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                data = metab)
summary(metab.m)

# extract model output and plot
metab.m2 <- data.frame(response = substr(row.names(metab.m$b), 9, 100), #remove the word response
                     estimate = metab.m$b, ci.lb = metab.m$ci.lb, ci.ub = metab.m$ci.ub)

metab.plot <- ggplot(metab.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#940845","#a9396a", "#bf6b8f", "#d49cb5", "#eaceda")) + # redish theme
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 6. Growth ##
grow <- subset(meta, trait == "Growth")

grow %>% # ES counts for growth response
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

grow %>% # n of study ID by growth response
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

grow.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = grow)
summary(grow.m)

# extract model output and plot
grow.m2 <- data.frame(response = substr(row.names(grow.m$b), 9, 100), #remove the word response
                      estimate = grow.m$b, ci.lb = grow.m$ci.lb, ci.ub = grow.m$ci.ub)

grow.plot <- ggplot(grow.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#941108","#a94139", "#bf706b", "#d4a09c")) + # orange theme
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 7. Immune function ##
immune <- subset(meta, trait == "Immune function")

immune %>% # ES counts for immune response
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))


immune %>% # n of study ID by immune response
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

immune.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                   random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                   control=list(optimizer="optim", optmethod="Nelder-Mead"), data = immune)
summary(immune.m)

# extract model output and plot
immune.m2 <- data.frame(response = substr(row.names(immune.m$b), 9, 100), #remove the word response
                        estimate = immune.m$b, ci.lb = immune.m$ci.lb, ci.ub = immune.m$ci.ub)

immune.plot <- ggplot(immune.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#945708", "#a97939", "#bf9a6b", "#d4bc9c", "#eaddce")) + # brown theme
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 8. Locomotion ##
loco <- subset(meta, trait == "Locomotion")

loco %>% # ES counts for locomotor response
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

loco %>% # n of study ID by locomotor response
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

# select 3 most well-sampled responses based on study ID
loco2 <- subset(loco, response %in% c("Ucrit","Speed","Muscle SERCA activity"))

loco.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = loco2)
summary(loco.m)

# extract model output and plot
loco.m2 <- data.frame(response = substr(row.names(loco.m$b), 9, 100), #remove the word response
                      estimate = loco.m$b, ci.lb = loco.m$ci.lb, ci.ub = loco.m$ci.ub)

loco.plot <- ggplot(loco.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#8B9408", "#a2a939", "#b9bf6b", "#d1d49c", "#e8eace")) + # lime theme
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 9. Plant-specific response ##
plant <- subset(meta, trait == "Plant-specific")

plant.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = plant)
summary(plant.m)

# extract model output and plot
plant.m2 <- data.frame(response = substr(row.names(plant.m$b), 9, 100), #remove the word response
                      estimate = plant.m$b, ci.lb = plant.m$ci.lb, ci.ub = plant.m$ci.ub)

locoplant.m <- rbind(loco.m2, plant.m2)
locoplant.m$trait <- c("Locomotion","Locomotion","Locomotion","Plant-specific","Plant-specific")

locoplant.plot <- ggplot(locoplant.m, aes(x = reorder(response, -estimate), y = estimate, colour = trait)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#8B9408","#459408")) + # lime and green theme
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 10. Reproduction ##
repro <- subset(meta, trait == "Reproduction")

repro %>% # ES counts for repro response
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

repro %>% # n of study ID by repro response
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

# select 5 most well-sampled responses based on study ID
repro2 <- subset(repro, response %in% c("Vtg level","Gonad mass","Offspring output", "E2 level", "Egg quantity"))

repro.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = repro2)
summary(repro.m)

# extract model output and plot
repro.m2 <- data.frame(response = substr(row.names(repro.m$b), 9, 100), #remove the word response
                      estimate = repro.m$b, ci.lb = repro.m$ci.lb, ci.ub = repro.m$ci.ub)

repro.plot <- ggplot(repro.m2, aes(x = reorder(response, -estimate), y = estimate, colour = response)) +
  mytheme + geom_hline(yintercept=0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#089411", "#39a941", "#6bbf70", "#9cd4a0", "#ceeacf")) + # light green theme
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

# 11. Thyroid response ##
thy <- subset(meta, trait == "Thyroid response")

thy %>% # ES counts for thyroid response
  group_by(response) %>% 
  summarise(count = n_distinct(effect_size_ID))

thy %>% # n of study ID bythyroid response
  group_by(response) %>% 
  summarise(count = n_distinct(study_ID))

thy.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                  random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                  data = thy)
summary(thy.m)

# extract model output and plot
thy.m2 <- data.frame(response = substr(row.names(thy.m$b), 9, 100), #remove the word response
                       estimate = thy.m$b, ci.lb = thy.m$ci.lb, ci.ub = thy.m$ci.ub)

# 12. Surivival response ##
surv <- subset(meta, trait == "Survival")

surv.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ response-1, 
                random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                data = surv)
summary(surv.m)

# extract model output and plot
surv.m2 <- data.frame(response = "Survival", #remove the word response
                     estimate = surv.m$b, ci.lb = surv.m$ci.lb, ci.ub = surv.m$ci.ub)

thy.m3 <- rbind(surv.m2, thy.m2)

thy.m3$trait <- c("Survival","Thyroid Response","Thyroid Response","Thyroid Response")

thy.plot <- ggplot(thy.m3, aes(x = reorder(response, -estimate), y = estimate, colour = trait)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 3.5, show.legend = FALSE) +
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  scale_color_manual(values = c("#089457", "#088B94")) + # green & turquoise theme
  xlab(NULL) + ylab(NULL) + 
  theme(axis.text.x = element_text(size = 8, angle = 10, hjust = 1))

library(cowplot)
plot_grid(abnorm.plot, behav.plot, cardio.plot, dev.plot, metab.plot, grow.plot,
          immune.plot, locoplant.plot, repro.plot, thy.plot,
          ncol = 2, align = "hv", axis = "tblr")

## Taxa overall ##--------------------------------------------------------------------
# within taxa alone
taxa.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ taxa-1, 
                  random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                  data = meta)

summary(taxa.m)

# Maginal R2 for taxa
fix <- var(as.numeric(as.vector(taxa.m$b) %*% t(as.matrix(taxa.m$X))))
R2m <- fix / (fix+sum(taxa.m$sigma2))
100 * R2m 

# I2
s2m <- sum(1/taxa.m$vi) * (taxa.m$k-1) / (sum(1/taxa.m$vi)^2 - sum((1/taxa.m$vi)^2)) 
s2t <- sum(taxa.m$sigma2) + s2m 

sum(taxa.m$sigma2) / s2t # Total heterogeneity
taxa.m$sigma2[1] / s2t # Species = non-phylogenetic variance
taxa.m$sigma2[2] / s2t # Group_ID = group (exp) variance
taxa.m$sigma2[3] / s2t # ES_ID = effect size level variance

# extract estimate and CIs
taxa.m3 <- data.frame(taxa = levels(meta$taxa),
                      estimate = taxa.m$b, ci.lb = taxa.m$ci.lb, ci.ub = taxa.m$ci.ub)
taxa.m3$taxa2 <- factor(taxa.m3$taxa, levels = c("Macroalgae","Phytoplankton","Angiosperms","Molluscs","Other invertebrates","Insects","Crustaceans","Amphibians","Fish"))
taxa.m3$domain <- c("Vertebrate","Autotroph","Invertebrate","Vertebrate","Invertebrate","Autotroph","Invertebrate","Invertebrate","Autotroph")

# plot in meta_phylo.R (ploteffect)

## Developmental period ##-------------------------------------
# within stage alone
levels(meta$development_stage)
meta$development_stage <- factor(meta$development_stage, 
                                 levels = c("Embryo", "Larvae", "Juvenile", "Adult"))

devel.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ development_stage-1, 
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = meta)
summary(devel.m)

# Maginal R2 for taxa
fix <- var(as.numeric(as.vector(devel.m$b) %*% t(as.matrix(devel.m$X))))
R2m <- fix / (fix+sum(devel.m$sigma2))
100 * R2m

# I2
s2m <- sum(1/devel.m$vi) * (devel.m$k-1) / (sum(1/devel.m$vi)^2 - sum((1/devel.m$vi)^2)) 
s2t <- sum(devel.m$sigma2) + s2m 

sum(devel.m$sigma2) / s2t # Total heterogeneity 
devel.m$sigma2[1] / s2t # Species = non-phylogenetic variance
devel.m$sigma2[2] / s2t # Group_ID = group (exp) variance 
devel.m$sigma2[3] / s2t # ES_ID = effect size level variance

# extract estimates and CIs, and plotting
devel.m3 <- data.frame(develop = levels(meta$development_stage),
                      estimate = devel.m$b, ci.lb = devel.m$ci.lb, ci.ub = devel.m$ci.ub)

ggplot(devel.m3, aes(x = reorder(develop, -estimate), y = estimate, colour = develop)) +
  mytheme + geom_hline(yintercept=0, linetype = "dashed") + 
  geom_point(size = 4, show.legend = FALSE) +
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub), size = 0.8, width=0.1, show.legend = FALSE) +
  xlab(NULL) + ylab("effect size (lnRR)") +
  scale_color_manual(values=c("#dadaeb","#9e9ac8", "#6a51a3", "#3f007d")) +
  theme(axis.text.x = element_text(size=8, angle = 20, hjust = 1))

# between taxa/development
crust.m2 <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~development_stage-1, 
                   random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                   data = meta, subset = taxa == "Crustaceans")
summary(amphi.m)

macro.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, 
                  random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                  data = meta, subset = taxa == "Phytoplankton")
levels(meta$taxa)
# extract estimate and CIs
macro.m2 <- data.frame(taxa = "Macroalgae", development = " ", 
                        estimate = macro.m$b, ci.lb = macro.m$ci.lb, ci.ub = macro.m$ci.ub)
phyto.m2 <- data.frame(taxa = "Phytoplankton", development = " ", 
                       estimate = phyto.m$b, ci.lb = phyto.m$ci.lb, ci.ub = phyto.m$ci.ub)
angio.m2 <- data.frame(taxa = "Angiosperms", development = substr(row.names(angio.m$b), 18, 100), 
                        estimate = angio.m$b, ci.lb = angio.m$ci.lb, ci.ub = angio.m$ci.ub)
mollusc.m2 <- data.frame(taxa = "Molluscs", development = substr(row.names(mollusc.m$b), 18, 100), #remove the word response
                         estimate = mollusc.m$b, ci.lb = mollusc.m$ci.lb, ci.ub = mollusc.m$ci.ub)
invert.m2 <- data.frame(taxa = "Other invertebrates", development = substr(row.names(invert.m$b), 18, 100),
                       estimate = invert.m$b, ci.lb = invert.m$ci.lb, ci.ub = invert.m$ci.ub)
insect.m2 <- data.frame(taxa = "Insects", development = substr(row.names(insect.m$b), 18, 100), 
                       estimate = insect.m$b, ci.lb = insect.m$ci.lb, ci.ub = insect.m$ci.ub)
crust.m2 <- data.frame(taxa = "Crustaceans", development = substr(row.names(crust.m$b), 18, 100), 
                        estimate = crust.m$b, ci.lb = crust.m$ci.lb, ci.ub = crust.m$ci.ub)
amphi.m2 <- data.frame(taxa = "Amphibians", development = substr(row.names(amphi.m$b), 18, 100), 
                        estimate = amphi.m$b, ci.lb = amphi.m$ci.lb, ci.ub = amphi.m$ci.ub)
fish.m2 <- data.frame(taxa = "Fish", development = substr(row.names(fish.m$b), 18, 100), 
                       estimate = fish.m$b, ci.lb = fish.m$ci.lb, ci.ub = fish.m$ci.ub)

devel.m4 <- rbind(macro.m2, phyto.m2, angio.m2, mollusc.m2, invert.m2,
                  insect.m2, crust.m2, amphi.m2, fish.m2)

devel.m4 <- rbind(angio.m2, mollusc.m2, invert.m2,
                  insect.m2, crust.m2, amphi.m2, fish.m2)

levels(devel.m4$development)
devel.m4$development <- factor(devel.m4$development, levels = c("Embryo", "Larvae", "Juvenile", "Adult", " "))

# plot in meta_phylo.R (plotdev)

## Between sex ##-----------------------------------------------------
sex <- meta %>%
  filter(sex == c("Male", "Female", "Both")) %>% 
  droplevels()

levels(sex$sex)

sex.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ sex-1,
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = sex)
summary(sex.m)

# Maginal R2 for taxa
fix < -var(as.numeric(as.vector(sex.m$b) %*% t(as.matrix(sex.m$X))))
R2m <- fix / (fix+sum(sex.m$sigma2))
100 * R2m

# I2
s2m <- sum(1/sex.m$vi) * (sex.m$k-1) / (sum(1/sex.m$vi)^2 - sum((1/sex.m$vi)^2)) 
s2t <- sum(sex.m$sigma2) + s2m 

sum(sex.m$sigma2) / s2t # Total heterogeneity
sex.m$sigma2[1] / s2t # Species = non-phylogenetic variance
sex.m$sigma2[2] / s2t # Group_ID = group (exp) variance
sex.m$sigma2[3] / s2t # ES_ID = effect size level variance

## Hierarchy ##--------------------------------------------------
hierh.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ hierarchy-1,
                random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                data = meta)
summary(hierh.m)

# Maginal R2 for level
fix <- var(as.numeric(as.vector(hierh.m$b) %*% t(as.matrix(hierh.m$X))))
R2m <- fix / (fix+sum(hierh.m$sigma2))
100 * R2m 

# I2
s2m <- sum(1/hierh.m$vi) * (hierh.m$k-1) / (sum(1/hierh.m$vi)^2 - sum((1/hierh.m$vi)^2)) 
s2t <- sum(hierh.m$sigma2) + s2m 

sum(hierh.m) / s2t # Total heterogeneity
hierh.m$sigma2[1] / s2t # Species = non-phylogenetic variance
hierh.m$sigma2[2] / s2t # Group_ID = group (exp) variance
hierh.m$sigma2[3] / s2t # ES_ID = effect size level variance

# plot by BPA concentration
# Reorder by food type
levels(meta$hierarchy)
meta$hierarchy <- factor(meta$hierarchy, levels = c("Organism","Organ","Tissue","Cellular","Macromolecules"))

BPA.plot <- ggplot(meta, aes(x = log10(corrected_dose), y = lnRR_direct, colour = hierarchy)) + 
  mytheme + geom_hline(yintercept=0, linetype = "dashed") +
  geom_point(cex = 3, alpha = 0.5, show.legend = FALSE) + 
  geom_smooth(method = lm, formula = y ~ poly(x, 2), se = FALSE, 
              size = 1.5, show.legend = FALSE) +
  scale_color_manual(values=c("#bd0026", "#f03b20", "#fd8d3c", "#feb24c", "#fed976")) +
  facet_grid(hierarchy ~ ., scales='free')

plot_grid(level.plot, BPA.plot, ncol = 2, 
          align = "hv", axis = "lrtb", labels=c("C", "D"))

## Additional variables temp, BPA conc ##---------------------------------------------
cont.m <- rma.mv(yi = lnRR_direct, V = lnRR_variance2, mod = ~ temperature + corrected_dose + time_day, 
                 random = list(~1|scientific_name, ~1|study_ID, ~1|effect_size_ID), 
                 data = meta)


summary(cont.m)

fix <- var(as.numeric(as.vector(cont.m$b) %*% t(as.matrix(cont.m$X))))
R2m <- fix / (fix + sum(cont.m$sigma2))
100 * R2m # account 23.14% of varience


## Acknowledgements ##----------------------------------------------------------------------
# R codes were obtained and modified from Farquharson (https://www.nature.com/articles/s41467-018-03500-9#Sec17), 
# and Nakagawa (https://osf.io/jwx6d/),
# thank you authors for making these codes freely avaliable online!
