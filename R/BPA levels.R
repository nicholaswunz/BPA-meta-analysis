# Author: Nicholas Wu (nicholas.wu.nz@gmail.com)
# Date: 23/04/2020
# R version: 3.5.1 -- "Feather Spray"
# Paper ID: BPA meta-analysis
# Description: Plotting BPA levels in taxa and the envrironment, produce global distribution map and analysis.

# install and load packages
install.packages("")
library(ggplot2) # use ggplot2 method of plotting
library(cowplot) # organise figures
library(dplyr) # for data manipulation like pipes %>%
#library(plyr) # for ddply summary stats
#detach("package:plyr", unload=TRUE) # clashes with dplyr at times
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
setwd('C:/Users/chwu6540/Dropbox (Personal)/Zebrafish physiology - Sydney/Meta analysis/BPA level data')

## Plot BPA level across taxa  ##--------------------------------------------------------------
# load data
organ <- read.csv("BPA level taxa.csv")
str(organ)
length(unique(organ$scientific_name))# number of species

organ %>%
  group_by(taxa) %>%
  summarise(count = n_distinct(scientific_name))

organ %>%
  group_by(taxa) %>%
  summarise(count = n_distinct(BPA))

# Compare different detection limit 
# at detection limit
den.plot1 <- ggplot(organ) + mytheme +
  geom_density(aes(log(BPA_zero)), colour = "grey", size = 1) + # replaced with zero
  geom_density(aes(log(estimated_BPA)), colour = "grey", size = 1) + # DL divided by 2
  geom_density(aes(log(estimated_BPA_2)), colour = "grey", size = 1) + # DL divided by sqrt of 2
  geom_density(aes(log(BPA)), colour = "#084594", size = 1.5) + # at detection limit
  ggtitle("At detection limit") + xlab(NULL)
  
# replaced with zero
den.plot2 <- ggplot(organ) + mytheme +
  geom_density(aes(log(BPA)), colour = "grey", size = 1) + # at detection limit
  geom_density(aes(log(estimated_BPA)), colour = "grey", size = 1) + # DL divided by 2
  geom_density(aes(log(estimated_BPA_2)), colour = "grey", size = 1) + # DL divided by sqrt of 2
  geom_density(aes(log(BPA_zero)), colour = "#94088B", size = 1.5) + # replaced with zero
  ggtitle("Detection limit replaced with zero") + xlab(NULL)

# DL divided by 2
den.plot3 <- ggplot(organ) + mytheme +
  geom_density(aes(log(BPA)), colour = "grey", size = 1) + # at detection limit
  geom_density(aes(log(BPA_zero)), colour = "grey", size = 1) + # replaced with zero
  geom_density(aes(log(estimated_BPA_2)), colour = "grey", size = 1) + # DL divided by sqrt of 2
  geom_density(aes(log(estimated_BPA)), colour = "#945708", size = 1.5) + # DL divided by 2
  ggtitle("Detection limit divided by two") + xlab(NULL)

# DL divided by sqrt of 2
den.plot4 <- ggplot(organ) + mytheme +
  geom_density(aes(log(BPA)), colour = "grey", size = 1) + # at detection limit
  geom_density(aes(log(BPA_zero)), colour = "grey", size = 1) + # replaced with zero
  geom_density(aes(log(estimated_BPA)), colour = "grey", size = 1) + # DL divided by 2
  geom_density(aes(log(estimated_BPA_2)), colour = "#089411", size = 1.5) + # DL divided by sqrt of 2
ggtitle("Detection limit divided by square root of two") + xlab(NULL)

plot_grid(den.plot1, den.plot2, den.plot3, den.plot4, ncol = 2, align = "vh", axis = "lr", labels = c("A", "B", "C", "D"))

organ$logBPA <- log(organ$BPA + 1)
organ$logBPA_zero <- log(organ$BPA_zero + 1)
organ$logestimated_BPA <- log(organ$estimated_BPA + 1)
organ$logestimated_BPA_2 <- log(organ$estimated_BPA_2 + 1)

install.packages("fitdistrplus")
library(fitdistrplus)

# calculate AIC for data distribution
fit.logBPA <- fitdist(organ$logBPA, "norm")
fit.logBPA$aic # 2760.513
plot(fit.logBPA)

fit.logBPA_zero <- fitdist(organ$logBPA_zero, "norm")
fit.logBPA_zero$aic # 2801.632
plot(fit.logBPA_zero)

fit.logestimated_BPA <- fitdist(organ$logestimated_BPA, "norm")
fit.logestimated_BPA$aic # 2760.513
plot(fit.logestimated_BPA)

fit.logestimated_BPA_2 <- fitdist(organ$logestimated_BPA_2, "norm")
fit.logestimated_BPA_2$aic # 2766.513
plot(fit.logestimated_BPA_2)

# combine all detection limit methods into one column
detect.lim <- melt(organ, id.vars = "paper_ID", measure.vars = c("logBPA","logBPA_zero","logestimated_BPA", "logestimated_BPA_2"))
head(detect.lim)

lim.m <- lm(value ~ variable, data = detect.lim)
anova(lim.m)

# calculate geometric mean BPA levels and SD
exp(mean(log(organ$BPA + 1)))
exp(sd(log(organ$BPA + 1)))

mean(organ$BPA)
sd(organ$BPA)

# plot taxa by BPA levels
fig1 <- ggplot(organ, aes(x = taxa, y = BPA, colour = taxa)) + mytheme +
  geom_jitter(position = position_dodge(0.8), cex = 3) +
  xlab(NULL) +
  ylab(expression("BPA concentrations (ng l"^"-1"*")")) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
  scale_y_sqrt() +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1","#6BAED6", "#2171B5", "#2171B5", "#08519C", "#08306B"))


ggplot(organ, aes(x = taxa, y = BPA, colour = taxa)) + mytheme +
  geom_flat_violin(alpha = 0.2, scale = "width") +
  geom_jitter(position = position_dodge(0.8), cex = 3) +
  ylab(expression("BPA concentrations (ng l"^"-1"*")")) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
  scale_y_log10() + xlab(NULL) +
  scale_color_manual(values = c("#C6DBEF", "#9ECAE1","#6BAED6", "#2171B5", "#2171B5", "#08519C", "#08306B"))


## Plot BPA levels in the enviroment ##--------------------------------------------------------------
#install.packages("data.table") # install package
library(data.table) # for freed function

# fread is an alternative to read.csv, but it brings things in as data.table objects
enviro <- fread("BPA level enviro.csv")
head(enviro)

# subset to surface waters only
enviro <- subset(enviro, specific_location == "surface water")
str(enviro)

# change to appropriate class
enviro$BPA <- as.numeric(enviro$BPA) # as numeric for BPA levels
enviro$BPA_zero <- as.numeric(enviro$BPA_zero) # as numeric for BPA levels
enviro$estimated_BPA <- as.numeric(enviro$estimated_BPA) # as numeric for BPA levels
enviro$estimated_BPA_2 <- as.numeric(enviro$estimated_BPA_2) # as numeric for BPA levels
enviro$environment <- as.factor(enviro$environment) # as factor

# Compare different detection limit 
# at detection limit
den.plot1 <- ggplot(enviro) + mytheme +
  geom_density(aes(log(BPA_zero)), colour = "grey", size = 1) + # replaced with zero
  geom_density(aes(log(estimated_BPA)), colour = "grey", size = 1) + # DL divided by 2
  geom_density(aes(log(estimated_BPA_2)), colour = "grey", size = 1) + # DL divided by sqrt of 2
  geom_density(aes(log(BPA)), colour = "#084594", size = 1.5) + # at detection limit
  ggtitle("At detection limit") + xlab(NULL)

# replaced with zero
den.plot2 <- ggplot(enviro) + mytheme +
  geom_density(aes(log(BPA)), colour = "grey", size = 1) + # at detection limit
  geom_density(aes(log(estimated_BPA)), colour = "grey", size = 1) + # DL divided by 2
  geom_density(aes(log(estimated_BPA_2)), colour = "grey", size = 1) + # DL divided by sqrt of 2
  geom_density(aes(log(BPA_zero)), colour = "#94088B", size = 1.5) + # replaced with zero
  ggtitle("Detection limit replaced with zero") + xlab(NULL)

# DL divided by 2
den.plot3 <- ggplot(enviro) + mytheme +
  geom_density(aes(log(BPA)), colour = "grey", size = 1) + # at detection limit
  geom_density(aes(log(BPA_zero)), colour = "grey", size = 1) + # replaced with zero
  geom_density(aes(log(estimated_BPA_2)), colour = "grey", size = 1) + # DL divided by sqrt of 2
  geom_density(aes(log(estimated_BPA)), colour = "#945708", size = 1.5) + # DL divided by 2
  ggtitle("Detection limit divided by two") + xlab(NULL)

# DL divided by sqrt of 2
den.plot4 <- ggplot(enviro) + mytheme +
  geom_density(aes(log(BPA)), colour = "grey", size = 1) + # at detection limit
  geom_density(aes(log(BPA_zero)), colour = "grey", size = 1) + # replaced with zero
  geom_density(aes(log(estimated_BPA)), colour = "grey", size = 1) + # DL divided by 2
  geom_density(aes(log(estimated_BPA_2)), colour = "#089411", size = 1.5) + # DL divided by sqrt of 2
  ggtitle("Detection limit divided by square root of two") + xlab(NULL)

plot_grid(den.plot1, den.plot2, den.plot3, den.plot4, ncol = 2, align = "vh", axis = "lr", labels = c("A", "B", "C", "D"))

enviro$logBPA <- log(enviro$BPA + 1)
enviro$logBPA_zero <- log(enviro$BPA_zero + 1)
enviro$logestimated_BPA <- log(enviro$estimated_BPA + 1)
enviro$logestimated_BPA_2 <- log(enviro$estimated_BPA_2 + 1)

# calculate AIC for data distribution
fit.logBPA <- fitdist(enviro$logBPA, "norm")
fit.logBPA$aic # AIC 5953.803
plot(fit.logBPA)

fit.logBPA_zero <- fitdist(enviro$logBPA_zero, "norm")
fit.logBPA_zero$aic # AIC 6062.099
plot(fit.logBPA)

fit.logestimated_BPA <- fitdist(enviro$logestimated_BPA, "norm")
fit.logestimated_BPA$aic # AIC 5975.526
plot(fit.logBPA)

fit.logestimated_BPA_2 <- fitdist(enviro$logestimated_BPA_2, "norm")
fit.logestimated_BPA_2$aic # AIC 5964.659
plot(fit.logBPA)

# combine all detection limit methods into one column
detect.lim <- melt(enviro, id.vars = "paper_ID", measure.vars = c("logBPA","logBPA_zero","logestimated_BPA", "logestimated_BPA_2"))
head(detect.lim)

lim.m <- lm(value ~ variable, data = detect.lim)
anova(lim.m)

# some basic data summary functions
length(unique(enviro$country))
range(enviro$BPA)
mean(enviro$BPA)
sd(enviro$BPA)

# geometric mean and SD
exp(mean(log(enviro$BPA + 1)))
exp(sd(log(enviro$BPA + 1)))


# graphing the geometric mean values per environment
enviromean <- enviro %>% 
  group_by(environment) %>% 
  summarise(mean = exp(mean(log(BPA + 1))),
            sd = exp(sd(log(BPA + 1))),
            max = max(BPA),
            min = min(BPA),
            n = length(BPA))

enviro$environment <- factor(enviro$environment, levels = c("freshwater","estuary","marine"))
# geom_flat-violin from: https://zenodo.org/record/1421371#.XoVahKgzYUF

 fig2 <- ggplot() + mytheme +
  geom_flat_violin(data = enviro, aes(x = environment, y = BPA),
                   alpha = 0.2, colour = "#e0e0e0") +
  geom_jitter(data = enviro, aes(x = environment, y = BPA), 
              position = position_dodge(0.8), cex = 2, colour = "#e0e0e0") +
  geom_pointrange(data = enviromean, aes(x = environment, y = mean, 
                                       ymin = mean-sd, ymax = mean+sd), 
                  size = 0.9, colour = "#084594") +
  ylab(expression("BPA concentrations (ng l"^"-1"*")")) +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
  scale_y_log10() + xlab(NULL)
 
# graphing the geometric mean values per country
countrymean <- enviro %>% 
  group_by(country) %>% 
  summarise(mean = exp(mean(log(BPA + 1))),
            sd = exp(sd(log(BPA + 1))),
            max = max(BPA),
            min = min(BPA + 1),
            n = length(BPA))

# reorder country list by decending geometric mean
country <- reorder(countrymean$country, -countrymean$mean)
countrymean$country <- factor(countrymean$country, levels = levels(country)) # for geomean data
levels(countrymean$country)

enviro$country <- factor(enviro$country, levels = levels(country)) # for raw data
levels(enviro$country)

fig3 <- ggplot(countrymean, aes(x = country, y = mean, colour = log10(mean))) +
  mytheme + 
  #geom_point(aes(y = max), colour = "grey", stat = "identity", size = 2) +
  #geom_linerange(aes(ymin = min, ymax = max), colour = "grey", linetype = "dotted") +
  geom_point(data = enviro, aes(y = BPA), colour = "#e0e0e0", size = 2) +
  geom_point(stat = "identity", size = 3) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.3) +
  theme(axis.text.x = element_text(size = 7, angle = 45, hjust = 1)) +
  scale_y_continuous(trans = 'log10') +
  xlab(NULL) + ylab(expression("Mean BPA concentrations (ng l"^"-1"*")")) +
  geom_text(aes(label = n), vjust=-3, size = 2) +
  scale_colour_gradient2(low = "#deebf7", mid = "#6baed6", high = "#084594") 


# combine fig1 and fig2
set1 <- plot_grid(fig1 + theme(legend.position="none"),
                  fig2 + theme(legend.position="none"), 
                  nrow = 1, align = "vh", axis = "lr", labels = c("A", "B"))

# combine fig1/2 and fig 3
plot_grid(set1, fig3 + theme(legend.position="none"),
          ncol = 1, align = "v", axis = "l", labels = c("A", "B"))


## BPA maps #---------------------------------------------------------------------
# load map
#install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel","ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
library(sf) # standardised way to encode spatial vectors
library(rnaturalearth) # load map
library(rnaturalearthdata) # load data for map

# create dataframe of countries
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# upload mean BPA levels per country
worldBPA <- read.csv("world_2.csv") # load BPA levels that match the world.csv (include NA)
head(worldBPA)
merged <- merge(world, worldBPA, by.x = "name", by.y = "name")

# plot global BPA levels (geometric)
map1 <- ggplot(data = merged) + mytheme +
  geom_sf(aes(fill = geo_BPA),colour = NA, size = 0.2) +
  geom_path(data = coastlines2, aes(x = long, y = lat, group = group), size = 0.2) +
  geom_path(data = rivers50, aes(x = long, y = lat, group = group), colour = "#6a7980", size = 0.2) +
  xlab(NULL) + ylab(NULL) + theme(panel.background = element_rect(fill = "#DEE5EB", colour = "black")) +
  scale_fill_distiller(palette = "Blues", direction = 1, trans = "log10", na.value = "#F1EEE7")

## Plastic pollution map ##-------------------------------------------------------
# Load packages into the workspace
library(rgdal) # Bindings for the 'Geospatial' Data Abstraction Library e.g. readOGR(), spTransform()
library(maptools) # for handling spatial objects
library(raster) # Geographic Data Analysis and Modeling
library(rgeos) # Interface to Geometry Engine
library(broom) # if you plot with ggplot and need to turn sp data into dataframes

# read plastic river shape file
plastic <- readOGR("PlasticRiverInput/PlasticRiverInputs.shp")
class(plastic) # check class type
extent(plastic)
crs(plastic) # coordinate system
plot(plastic)

# For coastline only shape file
# download the data
#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip",
#destfile = 'coastlines.zip')

# unzip the file
#unzip(zipfile = "coastlines.zip",
      #exdir = 'ne-coastlines-10m')

# from Lebreton et al. (2017) supplementary file
coastlines <- readOGR("PlasticRiverInput/ne-coastlines-10m/ne_10m_coastline.shp")
coastlines2 <- SpatialLinesDataFrame(coastlines, coastlines@data)

# check classes of all variables in spatial dataset
sapply(plastic@data, class)

# convert spatial data to dataframe
world2 <- data.frame(plastic)
str(world2)

# download rivers from R natural earth
rivers50 <- ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
sp::plot(rivers50)

# Plot inadequately-managed plastic wastev(2010, %) and average plastic input into rivers (tonne per year)
world2$i_mid1 <- world2$i_mid + 0.0000001 # for setting alpha range (doesnt work for NA's)

# plot global plastic waste
map2 <- ggplot() + mytheme +
  geom_point(data = world2, aes(x = coords.x1, y = coords.x2, colour = log10(i_mid), alpha = log10(i_mid1)), size = 5) +
  scale_alpha(range = c(0, 1), guide = FALSE) +
  geom_sf(data = merged, aes(fill = inadequat_waste), colour = NA, alpha = 0.9, size = 0.2) +
  geom_path(data = coastlines2, aes(x = long, y = lat, group = group), size = 0.2) +
  geom_path(data = rivers50, aes(x = long, y = lat, group = group), colour = "#806b6a", size = 0.2) +
  scale_colour_distiller(palette = "Oranges", direction = 1) +
  scale_fill_distiller(palette = "Oranges", direction = 1, na.value = "#F1EEE7") +
  xlab(NULL) + ylab(NULL) + theme(panel.background= element_rect(fill = "#DEE5EB", colour = "black"))

plot_grid(map1, map2, ncol = 1, align = "v", axis ="bt", labels=c("C", "D"))

# BPA plastic relationship
str(worldBPA)
worldBPA$inadequat_waste <- as.numeric(worldBPA$inadequat_waste)
worldBPA$plastic_gen <- as.numeric(worldBPA$plastic_gen)

library(ggpmisc) # for stat_poly_eq

ggplot(worldBPA, aes(x = log10(plastic_gen), y = log10(geo_BPA), colour = plastic_gen)) + 
  stat_smooth(method = 'lm', formula = 'y ~ x', colour = "#084594", fill = "#C6DBEF") +
  geom_point(size = 2) + mytheme +
  stat_poly_eq(formula = 'y ~ x', 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE) +
  scale_color_distiller(palette = "Blues", direction = 1, trans = "log10") +
  xlim(4,8) + geom_text(aes(label = name),hjust = -0.5, vjust = -0.5, size = 2) +
  theme(legend.position="none")
  

# lm for total plastic waste gen
plastic.md <- lm(log10(geo_BPA + 1) ~ log10(plastic_gen + 1), data = worldBPA)
summary(plastic.md)

# lm for mismanaged waste
waste.md <- lm(log10(geo_BPA + 1) ~ log10(inadequat_waste + 1), data = worldBPA)
summary(waste.md)
