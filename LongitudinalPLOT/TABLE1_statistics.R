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
library(rcompanion)

# PREPARE DATA ####
## Import data ####  
df_clean <- read.csv("tvals_cleansample_covars.csv") #todas las obs incluidas en el análisis

df_allIDs <- df_clean %>%
  distinct(ID, .keep_all = TRUE) #todos los sujetos incluidos en el análisis (solo 1 obs por sujeto)

df_long <- read.csv("subjects_finallist_DIT.csv") 
df_long <- df_long %>% drop_na(t)
df_long<- df_long %>% 
  semi_join(df_clean, by = "ID") #TABLA CON DUP y DIT pre sess1 

df_long <- df_long[!(df_long$ID == "1428" & df_long$session == 5), ]

df_baseline <- df_long %>% filter(t == 0) #datos de la ses1 de todos los sujetos incluidos en el análisis


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
  pairwise_fisher_test <- pairwiseNominalIndependence(contingency_table, fisher = TRUE, gtest = FALSE, chisq = FALSE, method = "bonferroni")
  print(pairwise_fisher_test)
}

rm(chi_square_test,male_percentage,pairwise_fisher_test, contingency_table)
# AGE ####
#1474 and 1509 only have ses 2 & 3 (no t = 0; missing APdose1 and PANSSTP_1)
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
}

rm(age_stats,kruskal_test,posthoc_test,shapiro_test_results)

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

rm(df_DUP, DUP_stats,wilcox_test)

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

rm(df_DIT, DIT_stats, wilcox_test)

# PANSS #####

df_PANSS <- df_long[df_long$Grupo!="HC",]
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

# APdose ####

#Fit the linear mixed model
model <- lmer(APdose_ ~ Resistance + (1 | ID), data = df_PANSS)
# Display the summary of the model
summary(model)
# Compute mean and standard deviation of PANSS_TP by group
APdose_stats <- df_PANSS %>%
  group_by(Grupo) %>%
  summarise(
    mean_APdose = format(round(mean(APdose_, na.rm = TRUE),2), nsmall = 2),
    sd_APdose = format(round(sd(APdose_, na.rm = TRUE),2), nsmall = 2)
  )
# Print the result
print(APdose_stats)

# IMPORT PANSS GEN and NEG ####
## Import new table ####  
sample_PANSS<-read.csv("subjects_May2024_PANSSTotal.csv")

#keep only columns and rows of interest
df_PANSSTOT_wide <- sample_PANSS[, c("ID", "Grupo", 
                 "PANSS_TP_1","PANSS_TP_2","PANSS_TP_3","PANSS_TP_4","PANSS_TP_5",
                 "PANSS_TN_1","PANSS_TN_2","PANSS_TN_3","PANSS_TN_4","PANSS_TN_5",
                 "PANSS_GEN_1","PANSS_GEN_2","PANSS_GEN_3","PANSS_GEN_4","PANSS_GEN_5",
                 "Excluir")]

#Factors: ID, Group (TRT, nonTRT), Sex
df_PANSSTOT_wide$ID <- as.factor(df_PANSSTOT_wide$ID)
df_PANSSTOT_wide$Grupo <- as.factor(df_PANSSTOT_wide$Grupo)
levels(df_PANSSTOT_wide$Grupo) 
levels(df_PANSSTOT_wide$Grupo) <- c("HC","nonTRS", "TRS")

#TO read numeric data:
df_PANSSTOT_wide$PANSS_TP_1  <- as.numeric(df_PANSSTOT_wide$PANSS_TP_1)
df_PANSSTOT_wide$PANSS_TP_2  <- as.numeric(df_PANSSTOT_wide$PANSS_TP_2)
df_PANSSTOT_wide$PANSS_TP_3  <- as.numeric(df_PANSSTOT_wide$PANSS_TP_3)
df_PANSSTOT_wide$PANSS_TP_4  <- as.numeric(df_PANSSTOT_wide$PANSS_TP_4)
df_PANSSTOT_wide$PANSS_TP_5  <- as.numeric(df_PANSSTOT_wide$PANSS_TP_5)

df_PANSSTOT_wide$PANSS_TN_1  <- as.numeric(df_PANSSTOT_wide$PANSS_TN_1)
df_PANSSTOT_wide$PANSS_TN_2  <- as.numeric(df_PANSSTOT_wide$PANSS_TN_2)
df_PANSSTOT_wide$PANSS_TN_3  <- as.numeric(df_PANSSTOT_wide$PANSS_TN_3)
df_PANSSTOT_wide$PANSS_TN_4  <- as.numeric(df_PANSSTOT_wide$PANSS_TN_4)
df_PANSSTOT_wide$PANSS_TN_5  <- as.numeric(df_PANSSTOT_wide$PANSS_TN_5)

df_PANSSTOT_wide$PANSS_GEN_1  <- as.numeric(df_PANSSTOT_wide$PANSS_GEN_1)
df_PANSSTOT_wide$PANSS_GEN_2  <- as.numeric(df_PANSSTOT_wide$PANSS_GEN_2)
df_PANSSTOT_wide$PANSS_GEN_3  <- as.numeric(df_PANSSTOT_wide$PANSS_GEN_3)
df_PANSSTOT_wide$PANSS_GEN_4  <- as.numeric(df_PANSSTOT_wide$PANSS_GEN_4)
df_PANSSTOT_wide$PANSS_GEN_5  <- as.numeric(df_PANSSTOT_wide$PANSS_GEN_5)
sapply(df_PANSSTOT_wide, class)

df_PANSSTOT_long <- pivot_longer(df_PANSSTOT_wide, 
                          cols =  starts_with("PANSS"),
                          names_to = c(".value", "session"),
                          names_pattern = "(\\D+)(\\d+)")
df_PANSSTOT_long$session <- as.integer(df_PANSSTOT_long$session)
df_PANSSTOT_long <- drop_na(df_PANSSTOT_long)

# Find rows in df_PANSSTOT that are not in df_PANSS
diff_PANSSTOT <- anti_join(df_PANSSTOT_long, df_PANSS, by = c("ID", "session"))
# Display the differences
print(diff_PANSSTOT)
# Remove rows in diff_PANSSTOT from df_PANSSTOT
df_PANSSTOT_cleaned <- anti_join(df_PANSSTOT_long, diff_PANSSTOT, by = c("ID", "session"))

df_PANSSTOT_cleaned$Resistance<- 0
df_PANSSTOT_cleaned$Resistance[df_PANSSTOT_cleaned$Grupo=="TRS"] <- 1
df_PANSSTOT_cleaned$Resistance <- as.factor(df_PANSSTOT_cleaned$Resistance)

# Fit the linear mixed model
model <- lmer(PANSS_TP_ ~ Resistance + (1 | ID), data = df_PANSSTOT_cleaned)
# Display the summary of the model
summary(model)
# Compute mean and standard deviation of PANSS_TP by group
PANSS_TP_stats <- df_PANSSTOT_cleaned %>%
  group_by(Grupo) %>%
  summarise(
    mean_PANSS_TP = format(round(mean(PANSS_TP_, na.rm = TRUE),2), nsmall = 2),
    sd_PANSS_TP = format(round(sd(PANSS_TP_, na.rm = TRUE),2), nsmall = 2)
  )
# Print the result
print(PANSS_TP_stats)


# Fit the linear mixed model
model <- lmer(PANSS_TN_ ~ Resistance + (1 | ID), data = df_PANSSTOT_cleaned)
# Display the summary of the model
summary(model)
# Compute mean and standard deviation of PANSS_TP by group
PANSS_TN_stats <- df_PANSSTOT_cleaned %>%
  group_by(Grupo) %>%
  summarise(
    mean_PANSS_TN = format(round(mean(PANSS_TN_, na.rm = TRUE),2), nsmall = 2),
    sd_PANSS_TN = format(round(sd(PANSS_TN_, na.rm = TRUE),2), nsmall = 2)
  )
# Print the result
print(PANSS_TN_stats)

# Fit the linear mixed model
model <- lmer(PANSS_GEN_ ~ Resistance + (1 | ID), data = df_PANSSTOT_cleaned)
# Display the summary of the model
summary(model)
# Compute mean and standard deviation of PANSS_TP by group
PANSS_GEN_stats <- df_PANSSTOT_cleaned %>%
  group_by(Grupo) %>%
  summarise(
    mean_PANSS_GEN = format(round(mean(PANSS_GEN_, na.rm = TRUE),2), nsmall = 2),
    sd_PANSS_GEN = format(round(sd(PANSS_GEN_, na.rm = TRUE),2), nsmall = 2)
  )
# Print the result
print(PANSS_GEN_stats)





library(patchwork)

# Assuming df_PANSSTOT_cleaned is already created as shown previously

# Create histogram for PANSS_TN_
p1 <- ggplot(df_PANSSTOT_cleaned, aes(x = PANSS_TN_, fill = Grupo)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
  labs(title = "PANSS_TN_ by Grupo", x = "PANSS_TN_", y = "Count", fill = "Grupo") +
  theme_minimal()

# Create histogram for PANSS_TP_
p2 <- ggplot(df_PANSSTOT_cleaned, aes(x = PANSS_TP_, fill = Grupo)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
  labs(title = "PANSS_TP_ by Grupo", x = "PANSS_TP_", y = "Count", fill = "Grupo") +
  theme_minimal()

# Create histogram for PANSS_GEN_
p3 <- ggplot(df_PANSSTOT_cleaned, aes(x = PANSS_GEN_, fill = Grupo)) +
  geom_histogram(binwidth = 1, position = "dodge", alpha = 0.7) +
  labs(title = "PANSS_GEN_ by Grupo", x = "PANSS_GEN_", y = "Count", fill = "Grupo") +
  theme_minimal()

# Combine the plots into a single column
combined_plot <- p1 / p2 / p3

# Display the combined plot
print(combined_plot)
