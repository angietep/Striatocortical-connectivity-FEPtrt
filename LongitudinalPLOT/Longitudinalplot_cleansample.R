# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
#setwd("/Users/brainsur/Desktop/GitHub_repos/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT")

library(tidyverse)
library(dplyr)
library(hrbrthemes)
#library(ggplot2)

# PREPARE DATA ####
## Import data ####  
df<-read.csv("tvals_cleansample_covars.csv")

# PLOT 1a: DURATION IN TREATMENT####

# Separate by groups (nonTRS and TRS)
nonTRS_data <- filter(df, group == "nonTRS")
TRS_data <- filter(df, group == "TRS")
HC_data <- filter(df, group == "HC")

# Then, reorder IDs based on the highest value of t_DIT for each subject
nonTRS_data <- nonTRS_data %>%
  arrange(desc(t_DIT)) %>%
  ungroup()

TRS_data <- TRS_data %>%
  arrange(desc(t_DIT)) %>%
  ungroup()

HC_data <- HC_data %>%
  arrange(desc(t_DIT)) %>%
  ungroup()


# Combine the separated and reordered data back together
organized_df <- rbind(TRS_data,nonTRS_data, HC_data)
organized_df$plotID<-as.numeric(row.names(organized_df))

# Convert ID to character and reorder based on plotID
organized_df$ID <- factor(organized_df$ID)
organized_df <- organized_df %>%
  group_by(group) %>%
  arrange(group, plotID) %>%
  mutate(ID = factor(ID, levels = rev(unique(ID)))) ## REMOVE REV TO HAVE THE INVERSE ORDER


## Plot ####

# Define the color palette manually or use another color palette function
color_palette <- c("#F8766D", "#00BFC4", "#7CAE00")  # Example colors, you can replace them with your preferred colors

ggplot(organized_df, aes(x = t_DIT)) +
  geom_line(aes(y = ID, group = ID), size = 0.5) +
  geom_point(aes(y = ID, color = as.factor(group)), size= 2) +
  labs(x = "Time between session (months)", y = "Participants", color = "Group") +
  scale_color_manual(name = "", 
                     values = color_palette,  # Change legend labels for color
                     labels = c("Healthy Controls", "Non Treatment-Resistant", "Treatment-Resistant")) +  # Set legend labels
  theme(
    legend.position = "top",  # Put legend inside figure (top-right corner)
    axis.text.y = element_blank(),  # Remove y-axis tick labels
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    legend.text = element_text(size = 14),  # Legend text size and style
    axis.title = element_text(size = 16, face = "bold"),  # Axis title size and style
    axis.text = element_text(size = 14),  # Axis tick labels size
    legend.title = element_text(size = 12, face = "bold")  # Legend title size and style
  )


df_baseline <- df[df$t_DIT==0,]
df_baseline$group <- factor(df_baseline$group, levels = c("nonTRS", "TRS","HC" ))
color_palette <- c("#00BFC4", "#7CAE00","#F8766D")  # Example colors, you can replace them with your preferred colors

ggplot(df_baseline, aes(x = age, fill = group)) +
  geom_histogram(position = "identity", alpha = 0.7, binwidth = 3) +  # Increase transparency for better visibility
  scale_fill_manual(values = color_palette, guide = "none") +  # Remove legend
  labs(x = "Age at first assessment", y = "Frequency") +
  facet_wrap(~ group, ncol = 1) +  # Display each group in a separate facet
  theme(
    axis.title = element_text(size = 16, face = "bold"),  # Axis title size and style
    axis.text = element_text(size = 14),  # Axis tick labels size
    strip.background = element_blank(),  # Remove facet labels background
    strip.text = element_blank()  # Remove facet labels
  )






                           

