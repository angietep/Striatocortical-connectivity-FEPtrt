rm(list = ls())

library(psych)
library(tidyverse)
library(forcats)

df <- read.csv("/Users/angeles/Documents/Research_data/Striatocortical_FEPtrs/cleansample_covars_R2_Med.csv")


#Total obs per group
table(df$group)

#Clozapina ####
table(df$Clozapina)
table(df$Clozapina,df$group)
# Convert to YES or NO 
df$Clozapina[!is.na(df$Clozapina)]<- 1
table(df$Clozapina, df$group)

df$group[df$Clozapina==1]<-"TRS_Clozapine"


#Olanzapina ####
df$Olanzapina[df$Olanzapina==""]<- NA
df$Olanzapina[df$Olanzapina=="Without doses"]<- "WOD"
table(df$Olanzapina, df$group) # 2 WOD
df$Olanzapina[df$Olanzapina=="WOD"]<- NA
# Convert to YES or NO 
df$Olanzapina[!is.na(df$Olanzapina)]<- 1
table(df$Olanzapina, df$group)

#Risperidona ####
table(df$Risperidona)
df$Risperidona[df$Risperidona==""]<- NA
df$Risperidona[df$Risperidona=="Without doses"]<-NA
table(df$Risperidona, df$group)
# Convert to YES or NO 
df$Risperidona[!is.na(df$Risperidona)]<- 1
table(df$Risperidona, df$group)

#Aripiprazol ####
table(df$Aripiprazol)
table(df$Aripiprazol, df$group)
# Convert to YES or NO 
df$Aripiprazol[!is.na(df$Aripiprazol)]<- 1
table(df$Aripiprazol, df$group)


#Quietapina ####
table(df$Quietapina)
df$Quietapina[df$Quietapina==""]<- NA
df$Quietapina[df$Quietapina=="Just SOS"]<- NA
table(df$Quietapina, df$group)
# Convert to YES or NO 
df$Quietapina[!is.na(df$Quietapina)]<- 1
table(df$Quietapina, df$group)


#Haloperidol  ####
table(df$Haloperidol)
table(df$Haloperidol,df$group)
# Convert to YES or NO 
df$Haloperidol[!is.na(df$Haloperidol)]<- 1
table(df$Haloperidol, df$group)

#Modecate ####
table(df$Modecate,df$group)
# Convert to YES or NO 
df$Modecate[!is.na(df$Modecate)]<- 1
table(df$Modecate, df$group)


#Amisulpride ####
table(df$Amisulpride, df$group)
# Convert to YES or NO 
df$Amisulpride[!is.na(df$Amisulpride)]<- 1
table(df$Amisulpride, df$group)

#Paliperidona ####
table(df$Paliperidona, df$group)
# Convert to YES or NO 
df$Paliperidona[!is.na(df$Paliperidona)]<- 1
table(df$Paliperidona, df$group)







# Define total N per group
group_totals <- c(NTR = 117, TR_beforeClozapine = 46, TR_afterClozapine = 15 )

# Define counts of people taking each medication by group
med_counts <- tribble(
  ~medication,    ~NTR, ~TR_beforeClozapine, ~TR_afterClozapine,
  "Olanzapina",       53,         28,              1,
  "Risperidona",      43,         27,              0,
  "Aripiprazol",      16,         10,              4,
  "Quetiapina",       5,          0,               0,
  "Haloperidol",      2,          1,               1,
  "Modecato",         0,          3,               0,
  "Amisulpride",      0,          1,               0,
  "Paliperidone",     2,          0,               0,
  "Clozapine",        0,          0,               15
)

med_counts <- tribble(
  ~medication,    ~NTR, ~TR_beforeClozapine, ~TR_afterClozapine,
  "Olanzapina",       53,         28,              1,
  "Risperidona",      43,         27,              0,
  "Aripiprazol",      16,         10,              4,
  "Quetiapina",       5,          0,               0,
  "Others",           4,          5,               1,
  "Clozapine",        0,          0,               15
)

med_counts <- tribble(
  ~medication,    ~NTR, ~TR_beforeClozapine, ~TR_afterClozapine,
  "Olanzapina",       53,         28,              1,
  "Risperidona",      43,         27,              0,
  "Aripiprazol",      16,         10,              4,
  "Others",           9,          5,               1,
  "Clozapine",        0,          0,               15
)


# Long format & percent calculation
df_plot <- med_counts %>%
  pivot_longer(cols = -medication, names_to = "group", values_to = "count") %>%
  mutate(total = group_totals[group],
         percentage = 100 * count / total) %>%
  group_by(medication) %>%
  mutate(avg_percentage = mean(percentage, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(medication = factor(medication,
                             levels = c("Clozapine", "Olanzapina", "Risperidona", 
                                        "Aripiprazol", "Others")))

# Plot
ggplot(df_plot, aes(x = group, y = percentage, fill = medication)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  coord_flip() +
  labs(title = "Percentage of Patients on Medication by Group",
       x = NULL,
       y = "Percentage (%)") +
  scale_fill_brewer(palette = "Set1", name = "Group") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

