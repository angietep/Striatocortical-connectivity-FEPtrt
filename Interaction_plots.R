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
sample<-read.csv("tvals_cleansample_covars.csv")

ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_timexTRS_cluster.1_size_721, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_timexTRS_cluster.2_size_536, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_time_cluster.1_size_644, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_SupVentralCaudate_time_cluster.2_size_504, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_VRPutamen_time_cluster.1_size_750, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_DCPutamen_time_cluster.1_size_441, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_InfVentralCaudate_time_cluster.1_size_595, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_DorsalCaudate_TRS_cluster.1_size_385, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_InfVentralCaudate_HC_cluster.1_size_1100, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()

ggplot(sample, aes(x = t_DIT, y = tvals_DorsalCaudate_HC_cluster.1_size_501, color = group, group = ID)) +
  geom_point() +
  geom_line(size = 0.5, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama()


