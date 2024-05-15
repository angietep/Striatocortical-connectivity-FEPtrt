# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
#setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
setwd("/Users/brainsur/Desktop/GitHub_repos/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT")

library(tidyverse)

# PREPARE DATA ####
## Import data ####  
sample<-read.csv("subjects_May2024.csv")
#keep only columns and rows of interest
df <- sample[, c("ID", "Grupo","Fecha_nacimiento","Sexo", 
                 "MRI_1","MRI_2","MRI_3","MRI_4", "MRI_5",
                 "APdose_1","APdose_2","APdose_3","APdose_4","APdose_5",
                 "PANSS_TP_1","PANSS_TP_2","PANSS_TP_3","PANSS_TP_4","PANSS_TP_5",
                 "DUP.months.","DIT_preSess1.days.",
                 "DIT_sess1toLastsess","Duration_of_illness.months.",
                 "Excluir")]

#keep ROWS up to subject 1517
#df <- df[1:117,] 
rm(sample) #remove original full df

## Change data classes ####
### PATIENTS ####
#Correct class and format of dates (Birth and MRI dates)
df$Fecha_nacimiento <- as.Date(df$Fecha_nacimiento, format = "%d/%m/%Y")
df$MRI_1 <- as.Date(df$MRI_1, format = "%d/%m/%Y")
df$MRI_2 <- as.Date(df$MRI_2, format = "%d/%m/%Y")
df$MRI_3 <- as.Date(df$MRI_3, format = "%d/%m/%Y")
df$MRI_4 <- as.Date(df$MRI_4, format = "%d/%m/%Y")
df$MRI_5 <- as.Date(df$MRI_5, format = "%d/%m/%Y")

#Factors: ID, Group (TRT, nonTRT), Sex
df$ID <- as.factor(df$ID)
df$Grupo <- as.factor(df$Grupo)
df$Sexo <- as.factor(df$Sexo)

# Rename the levels in Group: Resistente = TRS ; No_resistente = nonTRS
levels(df$Grupo) 
levels(df$Grupo) <- c("HC","nonTRS", "TRS")


#TO read numeric data:
# First) Replace commas with periods 
df$DUP.months. <- gsub(",", ".", df$DUP.months.)
df$DIT_preSess1.days. <- gsub(",", ".", df$DIT_preSess1.days.)
df$DIT_sess1toLastsess <- gsub(",", ".", df$DIT_sess1toLastsess)
df$Duration_of_illness.months. <- gsub(",", ".", df$Duration_of_illness.months.)
df$APdose_1 <- gsub(",", ".", df$APdose_1)
df$APdose_2 <- gsub(",", ".", df$APdose_2)
df$APdose_3 <- gsub(",", ".", df$APdose_3)
df$APdose_4 <- gsub(",", ".", df$APdose_4)
df$APdose_5 <- gsub(",", ".", df$APdose_5)
df$PANSS_TP_1 <- gsub(",", ".", df$PANSS_TP_1)
df$PANSS_TP_2 <- gsub(",", ".", df$PANSS_TP_2)
df$PANSS_TP_3 <- gsub(",", ".", df$PANSS_TP_3)
df$PANSS_TP_4 <- gsub(",", ".", df$PANSS_TP_4)
df$PANSS_TP_5 <- gsub(",", ".", df$PANSS_TP_5)

# Then) 
df$DUP.months. <- as.numeric(df$DUP.months.)
df$DIT_preSess1.days. <- as.numeric(df$DIT_preSess1.days.)
df$DIT_sess1toLastsess <- as.numeric(df$DIT_sess1toLastsess)
df$Duration_of_illness.months. <- as.numeric(df$Duration_of_illness.months.)
df$APdose_1 <- as.numeric(df$APdose_1)
df$APdose_2 <- as.numeric(df$APdose_2)
df$APdose_3 <- as.numeric(df$APdose_3)
df$APdose_4 <- as.numeric(df$APdose_4)
df$APdose_5 <- as.numeric(df$APdose_5)
df$PANSS_TP_1  <- as.numeric(df$PANSS_TP_1)
df$PANSS_TP_2  <- as.numeric(df$PANSS_TP_2)
df$PANSS_TP_3  <- as.numeric(df$PANSS_TP_3)
df$PANSS_TP_4  <- as.numeric(df$PANSS_TP_4)
df$PANSS_TP_5  <- as.numeric(df$PANSS_TP_5)

sapply(df, class)


## Exclude and remove NA ####
#PATIENTS
df <- df[!df$Excluir==1,]

df <- df[!is.na(df$Fecha_nacimiento), ] #excluded 0
df <- df[!is.na(df$Sexo), ] #excluded 0
df <- df[!is.na(df$Grupo), ] #excluded 2

#esta es la opción de incluirlos haciendo cero los NA
#df$DUP.months.[is.na(df$DUP.months.)] <- 0  


# PLOT 1a: DURATION IN TREATMENT####
## Define time-points for plot ####

# t0 = ses 1
df$t1 <- 0 
#time difference in "months" = weeks/4
df$t2<- df$t1 + as.numeric(difftime(df$MRI_2, df$MRI_1, units = 'weeks'))/4  
df$t3<- df$t1 + as.numeric(difftime(df$MRI_3, df$MRI_1, units = 'weeks'))/4
df$t4<- df$t1 + as.numeric(difftime(df$MRI_4, df$MRI_1, units = 'weeks'))/4
df$t5<- df$t1 + as.numeric(difftime(df$MRI_5, df$MRI_1, units = 'weeks'))/4


# Export the data frame to a CSV file
long_data <- pivot_longer(df, 
                          cols = starts_with("MRI") | starts_with("t") | starts_with("AP") | starts_with("PANSS"),
                          names_to = c(".value", "session"),
                          names_pattern = "(\\D+)(\\d+)")

long_data$PANSS_TP_[long_data$Grupo == "HC"] <- 0
long_data$APdose_[long_data$Grupo == "HC"] <- 0

long_data$Control[long_data$Grupo == "HC"] <- 1
long_data$Control[long_data$Grupo != "HC"] <- 0

long_data$Resistance[long_data$Grupo == "TRS"] <- 1
long_data$Resistance[long_data$Grupo != "TRS"] <- 0


write.csv(long_data, file = "subjects_finallist_DIT.csv", row.names = FALSE)






