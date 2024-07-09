# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
#setwd("/Users/brainsur/Desktop/striatconnTRT")
library(ggplot2)
library(geomtextpath)
library(ggsci)
library(dplyr)
library(tidyverse)
library(hrbrthemes)
library(viridis)

library(lme4)
library(emmeans)

# PREPARE DATA ####
## Import data ####  
sample<-read.csv("tvals_cleansample_covars.csv")

# Fit the model
#model_1 <- lmer(tvals_SupVentralCaudate_timexTRS_cluster.1_size_721 ~ age + sex + APdose + PANSS_TP + fdmean + HC + TRS + t_DIT + HC*t_DIT + TRS*t_DIT + (1 | ID), data = sample)
#model_2 <- lmer(tvals_SupVentralCaudate_timexTRS_cluster.2_size_536 ~ age + sex + APdose + PANSS_TP + fdmean + HC + TRS + t_DIT + HC*t_DIT + TRS*t_DIT + (1 | ID), data = sample)

# Fit the simpler model without interaction terms
model1_simpler <- lmer(tvals_SupVentralCaudate_timexTRS_cluster.1_size_721 ~ age + sex + APdose + PANSS_TP + fdmean + HC +  t_DIT + HC*t_DIT + (1 | ID), data = sample)
model2_simpler <- lmer(tvals_SupVentralCaudate_timexTRS_cluster.2_size_536 ~ age + sex + APdose + PANSS_TP + fdmean + HC + t_DIT + HC*t_DIT + (1 | ID), data = sample)

model3_simpler <- lmer(tvals_DorsalCaudate_TRS_cluster.1_size_385 ~ age + sex + APdose + PANSS_TP + fdmean + HC + t_DIT + HC*t_DIT + (1 | ID), data = sample)

# Extract residuals
sample$residuals1 <- resid(model1_simpler)
sample$residuals2 <- resid(model2_simpler)
sample$residuals3 <- resid(model3_simpler)

sample_PEP <- sample %>% filter(group!="HC")

# Plot using residuals
ggplot(sample_PEP, aes(x = t_DIT, y = residuals1, color = group, group = ID)) +
  geom_point(size = 3) +
  geom_line(size = 0.7, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama(labels = c("NTR", "TR")) +  # Change legend labels
  labs(
    x = "Time between sessions (months)",  # Replace with your desired x-axis label
    y = "Cluster's t-value"   # Replace with your desired y-axis label
  ) +
  theme(
    axis.title = element_text(size = 26, face = "bold"),  # Axis titles
    axis.text = element_text(size = 24, face = "bold"),   # Axis tick labels
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 24, face = "bold")  # Legend text
  )

ggplot(sample_PEP, aes(x = t_DIT, y = residuals2, color = group, group = ID)) +
  geom_point(size = 3) +
  geom_line(size = 0.7, alpha = 0.7,linetype = "dotted") +  # Thinner lines connecting points of the same subject
  geom_smooth(method = "lm", se = FALSE, size = 2, aes(group = group)) +  # Overall trend lines for each group
  theme_bw() +
  scale_color_jama(labels = c("NTR", "TR")) +  # Change legend labels
  labs(
    x = "Time between sessions (months)",  # Replace with your desired x-axis label
    y = "Cluster's t-value"   # Replace with your desired y-axis label
  ) +
  theme(
    axis.title = element_text(size = 26, face = "bold"),  # Axis titles
    axis.text = element_text(size = 24, face = "bold"),   # Axis tick labels
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 24, face = "bold")  # Legend text
  )

# Filter data for ses == 1
sample_ses1 <- sample_PEP %>% filter(ses == 1)

# Calculate means and CIs

summary_data <- sample_ses1 %>%
  group_by(group, t_DIT) %>%
  summarise(
    mean_residuals = mean(residuals3),
    ci_lower = mean_residuals - qt(0.975, df = n() - 1) * (sd(residuals3) / sqrt(n())),
    ci_upper = mean_residuals + qt(0.975, df = n() - 1) * (sd(residuals3) / sqrt(n()))
  )

summary_data$group[summary_data$group=="nonTRS"] <- "NTR"
summary_data$group[summary_data$group=="TRS"] <- "TR"

# Plot with CIs
plot_ci <- ggplot(summary_data, aes(x = group, y = mean_residuals, color = group)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), size= 1.5,width = 0.2) +
  theme_bw() +
  scale_color_jama(labels = c("NTR", "TR")) +  # Change legend labels
  labs(
    title = "Baseline (1st episode) Connectivity",
    x = "Groups",
    y = "Mean t-values (with CI)"
  ) +
  theme(
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 24, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 24, face = "bold"),
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5)  # Title styling
  )

# Display the plot
print(plot_ci)

