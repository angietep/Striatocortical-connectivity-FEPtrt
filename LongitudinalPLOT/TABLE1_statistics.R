# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
#setwd("/Users/brainsur/Desktop/GitHub_repos/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT")

library(tidyverse)
library(dplyr)
library(hrbrthemes)
library(rcompanion)
library(rstatix)
library(ggplot2)
library(car)
library(lme4)
library(lmerTest)
#library(ggplot2)

# PREPARE DATA ####
## Import data ####  
df_clean <- read.csv("cleansample_covars.csv") #todas las obs incluidas en el análisis
df_clean <- df_clean[df_clean$ID != 'C104', ] 
df_clean <- df_clean[df_clean$ID != 'C105', ] 

df_allIDs <- df_clean %>%
  distinct(ID, .keep_all = TRUE) #todos los sujetos incluidos en el análisis (solo 1 obs por sujeto)

df_long <- read.csv("subjects_finallist_DIT.csv") 
df_long <- df_long %>% drop_na(t)

df_long_filtered <- df_long %>% 
  semi_join(df_clean, by = "ID")

df_long_filtered <- df_long_filtered[!(df_long_filtered$ID == "1428" & df_long_filtered$session == 5), ]

df_baseline <- df_long_filtered %>% filter(t == 0) #datos de la ses1 de todos los sujetos incluidos en el análisis


# SEX ####
# Calculate the percentage of male in each group
male_percentage <- df_allIDs %>%
       group_by(group) %>%
       summarise(total = n(),
       males = sum(sex == "Masculino"),
       percentage_male = (males / total) * 100)

print(male_percentage)

# Create a contingency table for the chi-square test
contingency_table <- table(df_allIDs$group, df_allIDs$sex)

# Perform the chi-square test
chi_square_test <- chisq.test(contingency_table)
print(chi_square_test)

# Post-hoc pairwise comparisons if chi-square is significant
if (chi_square_test$p.value < 0.05) {
  # Using pairwise Fisher's exact test for post-hoc analysis
  library(rcompanion)
  pairwise_fisher_test <- pairwiseNominalIndependence(contingency_table, fisher = TRUE, gtest = FALSE, chisq = FALSE, method = "bonferroni")
  print(pairwise_fisher_test)
}

# AGE ####
#1474 and 1509 only have ses 2 & 3 (no t = 0; missing APdose1)
# find sess != 1 
# find those IDs in subjects_finallist and compute age as MRI - Fecha_NAC
# then mean and sd of age

df_allIDs$age[df_allIDs$ID=="1474"] <-as.integer(
  (as.Date(df_long$MRI_[df_long$ID == "1474" & df_long$t == 0]) - as.Date(df_long$Fecha_nacimiento[df_long$ID == "1474" & df_long$t == 0])) / 365
)
df_allIDs$age[df_allIDs$ID=="1509"] <-as.integer(
  (as.Date(df_long$MRI_[df_long$ID == "1509" & df_long$t == 0]) - as.Date(df_long$Fecha_nacimiento[df_long$ID == "1509" & df_long$t == 0])) / 365
)

# Calculate mean and standard deviation of age for each group
age_stats <- df_allIDs %>%
  group_by(group) %>%
  summarise(mean_age = mean(age, na.rm = TRUE),
            sd_age = sd(age, na.rm = TRUE))

print(age_stats)

shapiro_test_results <- df_allIDs %>%
  group_by(group) %>%
  summarise(p_value = shapiro.test(age)$p.value)

print(shapiro_test_results)
# Determine if any group deviates from normality
if (any(shapiro_test_results$p_value < 0.05)) {
  # Perform Kruskal-Wallis test if normality is violated
  kruskal_test <- kruskal_test(age ~ group, data = df_allIDs)
  print(kruskal_test)
  
  # Post-hoc Dunn test if Kruskal-Wallis is significant
  if (kruskal_test$p < 0.05) {
    posthoc_test <- dunn_test(age ~ group, data = df_allIDs, p.adjust.method = "bonferroni")
    print(posthoc_test)
  }
} else {
  # Perform ANOVA if normality is not violated
  anova_test <- aov(age ~ group, data = df_allIDs)
  anova_summary <- summary(anova_test)
  print(anova_summary)
  
  # Post-hoc Tukey test if ANOVA is significant
  if (anova_summary[[1]][["Pr(>F)"]][1] < 0.05) {
    posthoc_test <- TukeyHSD(anova_test)
    print(posthoc_test)
  }
}

df_allIDs <- df_clean %>%
  distinct(ID, .keep_all = TRUE) #todos los sujetos incluidos en el análisis (solo 1 obs por sujeto)


# DUP ####

df_DUP <- df_baseline %>% drop_na(DUP.months.)

DUP_stats <- df_DUP %>%
  group_by(Grupo) %>%
  summarise(mediand_DUP = format(round(median(DUP.months., na.rm = TRUE), 2), nsmall = 2),
            Q1_DUP = format(round(quantile(DUP.months.,0.25, na.rm = TRUE), 2), nsmall = 2),
            Q3_DUP = format(round(quantile(DUP.months.,0.75, na.rm = TRUE), 2), nsmall = 2))

print(DUP_stats)

df_DUP$Grupo <- as.factor(df_DUP$Grupo)
# Perform Wilcoxon rank-sum test
wilcox_test <- wilcox.test(DUP.months. ~ Grupo, data = df_DUP)
# Print the test results
print(wilcox_test)

rm(df_DUP)

# DIT pre sess1 ####

df_DIT <- df_baseline %>% drop_na(DIT_preSess1.days.)

DIT_stats <- df_DIT %>%
  group_by(Grupo) %>%
  summarise(
    mediand_DIT = format(round(median(DIT_preSess1.days., na.rm = TRUE), 2), nsmall = 2),
    Q1_DIT = format(round(quantile(DIT_preSess1.days., 0.25, na.rm = TRUE), 2), nsmall = 2),
    Q3_DIT = format(round(quantile(DIT_preSess1.days., 0.75, na.rm = TRUE), 2), nsmall = 2)
  )

print(DIT_stats)
sum(df_DIT$Grupo=="nonTRS")
sum(df_DIT$Grupo=="TRS")

df_DIT$Grupo <- as.factor(df_DIT$Grupo)
# Perform Wilcoxon rank-sum test
wilcox_test <- wilcox.test(DIT_preSess1.days. ~ Grupo, data = df_DIT)
# Print the test results
print(wilcox_test)

rm(df_DIT)

# PANSS #####

df_PANSS <- df_long_filtered[df_long_filtered$Grupo!="HC",]
df_PANSS <- df_PANSS[!(df_PANSS$ID == "1474" & df_PANSS$session == 1),] 
df_PANSS <- df_PANSS[!(df_PANSS$ID == "1509" & df_PANSS$session == 1),] 
df_PANSS$Resistance <- as.factor(df_PANSS$Resistance)

# Fit the linear mixed model
model <- lmer(PANSS_TP_ ~ Resistance + (1 | ID), data = df_PANSS)

# Display the summary of the model
summary(model)

# Compute mean and standard deviation of PANSS_TP by group
PANSS_TP_stats <- df_PANSS %>%
  group_by(Grupo) %>%
  summarise(
    mean_PANSS_TP = format(round(mean(PANSS_TP_, na.rm = TRUE),2), nsmall = 2),
    sd_PANSS_TP = format(round(sd(PANSS_TP_, na.rm = TRUE),2), nsmall = 2)
  )

# Print the result
print(PANSS_TP_stats)

