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
df<-read.csv("cleansample_covars.csv")
df <- df[df$ID != 'C104', ] 
df <- df[df$ID != 'C105', ] 


# LEER DATOS DE LA OTRA TABLA (SUBJECTS_DIT)
# Y CALCULAR EDADES Y SEXO INCLUYENDO A 1474 Y 1509 QUE NO TIENE APDOSE PARA SES 1


df_fup <- df %>%
     group_by(subject_id, group) %>%
     summarise(has_follow_up = n() > 1) %>%
     ungroup()

contingency_table <- table(df_fup$group, df_fup$has_follow_up)
print(contingency_table)

chi_square_test <- chisq.test(contingency_table)
print(chi_square_test)


contingency_table <- table(df$group[df$ses==1], df$sex[df$ses==1])
print(contingency_table)
chi_square_test <- chisq.test(contingency_table)


shapiro_test_HC <- shapiro.test(df$age[df$group == "HC" & df$ses==1])
shapiro_test_TRS <- shapiro.test(df$age[df$group == "TRS"  & df$ses==1])
shapiro_test_nonTRS <- shapiro.test(df$age[df$group == "nonTRS"  & df$ses==1])
print (shapiro_test_HC)


kruskal_test <- kruskal.test(age ~ group, data = df[df$ses==1, ])
print(kruskal_test)
pairwise_comparisons <- pairwise.wilcox.test(baseline_df$age, baseline_df$group, p.adjust.method = "bonferroni")
print(pairwise_comparisons)

