# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
#setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
setwd("/Users/brainsur/Desktop/striatconnTRT")
library(ggplot2)
library(geomtextpath)
library(ggsci)
# PREPARE DATA ####
## Import data ####  
sample<-read.csv("cleansample_covars.csv")


ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_interactionDITxgroup_cluster.3_size_839, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_interactionDITxgroup_cluster.1_size_2317, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.5) +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() 

ggplot(sample, aes(x = t_DIT, y = tvals_InfVentralCaudate_interactionDITxgroup_cluster.1_size_409, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.5) +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() 

ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_interactionDITxgroup_cluster.2_size_826, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.5) +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() 