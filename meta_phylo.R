# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 11/10/2019
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: BPA meta-analysis
# Description: Build phylogeny and estimate phylogenetic relatedness among species included in the meta analysis

# install and load packages
install.packages("")
library(rotl) # Interface to the 'Open Tree of Life' API.
library(ape) # deal with phylogenetic data

## Import dataset ##---------------------------------------------------------------------------------
#meta <- sublethal # load sublethal data from the meta_analysis.R 
head(meta) # check column headers
str(meta) # check strings
meta$scientific_name_OTL <- as.character(meta$scientific_name_OTL) # sort taxa names as character

# generating list of species
species <- sort(unique(meta$scientific_name_OTL))

## Formatting species data ##------------------------------------------------------------------------
# obtaining dataframe listing the Open Tree identifiers potentially matching our list of species.
taxa <- tnrs_match_names(names = species)
taxa # one species showed up "incertae sedis" "Cochlodinium polykrikoides"
# **rerun tnrs_match_names again**

# check if species list matcg OT identifier
taxa[taxa$approximate_match == TRUE,] # none so far

# retrieving phylogenetic relationships among taxa in the form of a trimmed sub-tree
tree <- tol_induced_subtree(ott_ids = ott_id(taxa),label_format = "name")

plot(tree, cex = 0.6, label.offset = 0.1, no.margin = TRUE) # plot tree

## Dealing with polytomies ##--------------------------------------------------------------------------
# If polytomies exist, the output will be `FALSE`, and vice versa.
is.binary.tree(tree)

# Use a randomization approach to deal with polytomies
set.seed(111) # making it replicable, at least for this version of R (i.e. v.3.5.1)
tree_random <- multi2di(tree, random=TRUE)

is.binary.tree(tree_random) # recheck output

plot(tree_random, cex = 0.6, label.offset = 0.1, no.margin = TRUE)

# exploring whether our tree covers all the species we wanted 
# it to include, and making sure that the species names in our 
# database match those in the tree. We use the following code.
tree_random$tip.label <- gsub("_"," ", tree_random$tip.label)
intersect(as.character(tree_random$tip.label), as.character(species))
setdiff(species, as.character(tree_random$tip.label)) # listed in the database but not in the tree
setdiff(as.character(tree_random$tip.label), species) # listed in the tree but not in the database

## Computing branch lengths ##----------------------------------------------------------------
# before we need to make sure that tree labels and database use the same nomenclature
setdiff(meta$scientific_name_OTL, as.character(tree_random$tip.label))
setdiff(as.character(tree_random$tip.label), meta$scientific_name_OTL)
tree_random.fixed <- tree_random

# gsub("in tree or old", "in database or new", string)
tree_random.fixed$tip.label <- gsub("Bufotes viridis","Bufo viridis", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Chlorolobion braunii","Monoraphidium braunii", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Navicula salinicola","Navicula incerta", tree_random.fixed$tip.label)
tree_random.fixed$tip.label <- gsub("Prorocentrum minimum (species in subkingdom SAR)","Prorocentrum minimum", tree_random.fixed$tip.label)

setdiff(meta$scientific_name_OTL, as.character(tree_random.fixed$tip.label))
setdiff(as.character(tree_random.fixed$tip.label),meta$scientific_name_OTL)
# all good!

# compute branch lengths of tree
phylo_branch <- compute.brlen(tree_random.fixed, method = "Grafen", power = 1)

# check tree is ultrametric
is.ultrametric(phylo_branch) # TRUE

## Phylogenetic matrix ##-------------------------------------------------------------------------
# matrix to be included in the models
phylo_cor <- vcv(phylo_branch, cor = T)
plot(phylo_branch)

# save matrix for future analyses
save(phylo_cor, file = "folder to save in")

## Plot phylogenetic tree ##------------------------------------------------------------------
# to install ggtree
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")
library(ggplot2)
library(ggtree) # plotting trees in ggplot2 format
library(dplyr)
library(cowplot) # organise figures

mytheme <- theme_bw() + {theme(panel.border = element_blank(), # Border around plotting area.
                              panel.grid.major = element_blank(), # Major grid lines blank
                              panel.grid.minor = element_blank(), # Minor grid lines blank
                              axis.line = element_line(colour = "black", size = 0.8), # axis line size
                              axis.ticks = element_line(colour = "black", size = 0.8),
                              axis.text = element_text(size = 10, colour = "black"), # axis text size
                              axis.title=element_text(size = 10), #axis title size
                              panel.background = element_rect(fill = "transparent"), # bg of the panel
                              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
                              legend.background = element_rect(fill = "transparent"), # get rid of legend bg
                              legend.box.background = element_rect(fill = "transparent", color = NA)) # get rid of legend panel bg)
}

# visualise nodes to group taxa
ggtree(phylo_branch) + geom_text(aes(label=node), hjust=-0.2, size = 2.5) + 
  geom_tiplab(size=3, hjust = -0.2) + xlim_tree(1.5)

# group by node numbers as domain
cls <- list(Vertebrate = c(48:82),
            Invertebrate = c(20:47, 83:87),
            Autotroph = c(1:19, 88))

# group by node numbers as taxa
taxa <- list(amphibian = c(51:60),
             angiosperm = c(1:3),
             crustacean = c(38:50),
             fish = c(61:85),
             insect = c(36:37),
             macroalgae = c(15,4),
             mollusc = c(24:32),
             "other invertebrates" = c(86:90,23, 33:35),
             phytoplankton = c(5:14, 16:22, 91))

# merge group values into main tree
tree2 <- groupOTU(phylo_branch, cls)

ggtree(tree2, aes(color = group), size = 1) + geom_tiplab(size = 3) +
  scale_color_manual(values = c("#9ecae1", "#4292c6", "#084594")) +
  xlim_tree(1.5)

# plot tree by domain
plottree <- ggtree(tree2, aes(color = group), size = 1)+ xlim_tree(1.2) +
  geom_cladelabel(node = 48:82, label = "Vertebrates", align = T, barsize = 1, fontsize = 3, vjust = -0.5) +
  geom_cladelabel(node = c(20:47, 83:87), label = "Invertebrates", align = T, barsize = 1, fontsize = 3, vjust = -0.5) +
  geom_cladelabel(node = c(1:19, 88), label = "Autotrophs", align = T, barsize = 1, fontsize = 3, vjust = -0.5) +
  scale_color_manual(values = c("#9ecae1", "#4292c6", "#084594")) +
  coord_flip()


## plot with effects size (grouped taxa) ##-------------------------------------------------------------------

# function to sqrt negative values
library(scales)
S_sqrt <- function(x){sign(x)*sqrt(abs(x))}
IS_sqrt <- function(x){x^2*sign(x)}
S_sqrt_trans <- function() trans_new("S_sqrt",S_sqrt,IS_sqrt)

# plot mean taxa effect size
ploteffect <- ggplot(taxa.m3, aes(x = taxa2, y = estimate, colour = domain)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 4, show.legend = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, show.legend = FALSE) +
  xlab(NULL) + ylab("effect size (lnRR)") +
  scale_color_manual(values = c("#9ecae1", "#4292c6", "#084594")) +
  theme(axis.text.x = element_text(size = 8, angle = 20, hjust = 1))
  
# plot mean taxa effect size between development stage
plotdev <- ggplot(devel.m4, aes(x = taxa, y = estimate, colour = development)) +
  mytheme + geom_hline(yintercept = 0, linetype = "dashed") + 
  geom_point(size = 4, position = position_dodge(width = 0.8), show.legend = FALSE)+
  geom_errorbar(aes(ymin = ci.lb, ymax = ci.ub), size = 0.8, width = 0.1, 
                position = position_dodge(width = 0.8), show.legend = FALSE) +
  xlab(NULL) + ylab("effect size (lnRR)") +
  scale_color_manual(values = c("#dadaeb", "#9e9ac8", "#6a51a3", "#3f007d", "black")) +
  theme(axis.text.x = element_text(size = 8, angle = 20, hjust = 1))

# development only small graph in meta_analysis.R

# combine effect size with phylogeny tree
plot_grid(ploteffect, plottree, plotdev, ncol = 1, 
          align = "v", axis = "rl", 
          rel_heights = c(1,0.7,1), labels = c("A", " ", "B"))


## Acknowledgements ##----------------------------------------------------------------------
# R code were obtained and  modified from ASanche-Tojar (https://github.com/ASanchez-Tojar/meta-analysis_of_variance), thank you authors!
