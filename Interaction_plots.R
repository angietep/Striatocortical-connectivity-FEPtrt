# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
#setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
setwd("/Users/brainsur/Desktop/striatconnTRT")
library(ggplot2)
library(geomtextpath)
library(ggsci)
library(dplyr)
# PREPARE DATA ####
## Import data ####  
sample<-read.csv("tvals_cleansample_covars.csv")

## Interaction PLOTS ####
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







## Demographics and statistics ####

df_last_session <- sample %>%
  group_by(ID) %>%
  filter(ses == max(ses)) %>%
  ungroup()

follow_ups <- df_last_session %>%
   filter(ses!=1) %>%
   ungroup()

boxplot(t_DIT ~ group, data = follow_ups, main = "Follow-ups", ylab = "Group", xlab = "Follow-up time (months)", col = "lightblue", horizontal = TRUE)

median(follow_ups$t_DIT[follow_ups$group == "HC"])
IQR(follow_ups$t_DIT[follow_ups$group == "HC"])
quantile(follow_ups$t_DIT[follow_ups$group == "HC"], 0.25)
quantile(follow_ups$t_DIT[follow_ups$group == "HC"], 0.75)

median(follow_ups$t_DIT[follow_ups$group == "TRS"])
IQR(follow_ups$t_DIT[follow_ups$group == "TRS"])
quantile(follow_ups$t_DIT[follow_ups$group == "TRS"], 0.25)
quantile(follow_ups$t_DIT[follow_ups$group == "TRS"], 0.75)

median(follow_ups$t_DIT[follow_ups$group == "nonTRS"])
IQR(follow_ups$t_DIT[follow_ups$group == "nonTRS"])
quantile(follow_ups$t_DIT[follow_ups$group == "nonTRS"], 0.25)
quantile(follow_ups$t_DIT[follow_ups$group == "nonTRS"], 0.75)
