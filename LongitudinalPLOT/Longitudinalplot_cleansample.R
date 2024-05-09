# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
#setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
setwd("/Users/brainsur/Desktop/GitHub_repos/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT")

library(tidyverse)
library(dplyr)

# PREPARE DATA ####
## Import data ####  
df<-read.csv("cleansample_covars.csv")

# PLOT 1a: DURATION IN TREATMENT####

# Separate by groups (nonTRT and TRT)
nonTRT_data <- filter(df, group == "nonTRT")
TRT_data <- filter(df, group == "TRT")

# Then, reorder IDs based on the highest value of t_DIT for each subject
nonTRT_data <- nonTRT_data %>%
  arrange(desc(t_DIT)) %>%
  ungroup()

TRT_data <- TRT_data %>%
  arrange(desc(t_DIT)) %>%
  ungroup()

# Combine the separated and reordered data back together
organized_df <- rbind(nonTRT_data, TRT_data)

organized_df$plotID<-as.numeric(row.names(organized_df))

# Convert ID to character and reorder based on plotID
organized_df$ID <- factor(organized_df$ID)
organized_df <- organized_df %>%
  group_by(group) %>%
  arrange(group, plotID) %>%
  mutate(ID = factor(ID, levels = rev(unique(ID)))) ## REMOVE REV TO HAVE THE INVERSE ORDER

specific_ticks <- c(12, 12*2, 12*4, 12*6, 12*7)
specific_labels <- c("","","","","")
## Plot ####
ggplot(organized_df, aes(x = t_DIT)) +
  geom_line(aes(y=ID, group = ID ), size = .3) +
  geom_point(aes(y=ID, color = as.factor(group),shape = factor(ses)), size=2) +
  labs(x = "Duration of treatment (months)", y = "Subject ID", color = "Group", shape= "Session") +
  theme(legend.position = "top", axis.text.y = element_text(size = 5)) +
  geom_density(aes(y = -..density..*1000), fill= "#69b3a2") +
  geom_label(aes(x=75, y=-10, label="density")) + 
  scale_x_continuous(breaks = c(specific_ticks, 0, 25, 50, 75, 100),
                     labels = c(specific_labels, 0, 25, 50, 75, 100)) +
  geom_text(aes(x=12, y=-20, label="1y"), size=3) +
  geom_text(aes(x=12*2, y=-20, label="2y"), size=3) +
  geom_text(aes(x=12*4, y=-20, label="4y"), size=3) +
  geom_text(aes(x=12*6, y=-20, label="6y"), size=3) +
  geom_text(aes(x=12*7, y=-20, label="7y"), size=3)



