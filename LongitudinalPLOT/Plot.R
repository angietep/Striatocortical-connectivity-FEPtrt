# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list= ls())# ctrl + L to clear console
setwd("~/Documents/GitHub/Striatocortical-connectivity-FEPtrt/LongitudinalPLOT/")
library(tidyverse)

# PREPARE DATA ####
## Import data ####  
sample<-read.csv("subjects_PEP_Jan2024.csv")
#keep only columns and rows of interest
df <- sample[, c("ID", "Grupo","Fecha_nacimiento","Sexo", 
                 "MRI_1","MRI_2","MRI_3","MRI_4",
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

#Factors: ID, Group (TRT, nonTRT), Sex
df$ID <- as.factor(df$ID)
df$Grupo <- as.factor(df$Grupo)
df$Sexo <- as.factor(df$Sexo)

# Rename the levels in Group: Resistente = TRT ; No_resistente = nonTRT
levels(df$Grupo) 
levels(df$Grupo) <- c("nonTRT", "TRT")


#TO read numeric data:
# First) Replace commas with periods 
df$DUP.months. <- gsub(",", ".", df$DUP.months.)
df$DIT_preSess1.days. <- gsub(",", ".", df$DIT_preSess1.days.)
df$DIT_sess1toLastsess <- gsub(",", ".", df$DIT_sess1toLastsess)
df$Duration_of_illness.months. <- gsub(",", ".", df$Duration_of_illness.months.)
# Then) 
df$DUP.months. <- as.numeric(df$DUP.months.)
df$DIT_preSess1.days. <- as.numeric(df$DIT_preSess1.days.)
df$DIT_sess1toLastsess <- as.numeric(df$DIT_sess1toLastsess)
df$Duration_of_illness.months. <- as.numeric(df$Duration_of_illness.months.)

sapply(df, class)


## Exclude and remove NA ####
#PATIENTS
df <- df[!df$Excluir==1,]

df <- df[!is.na(df$Fecha_nacimiento), ]
df <- df[!is.na(df$Sexo), ]
df <- df[!is.na(df$Grupo), ]
df <- df[!is.na(df$DIT_preSess1.days.), ]


#esta es la opción de incluirlos haciendo cero los NA
#df$DUP.months.[is.na(df$DUP.months.)] <- 0  


# PLOT 1a: DURATION IN TREATMENT####
## Define time-points for plot ####

# t0 = 0 : inicio de enfermedad
# t1: time from illness start to first MRI scan
#t1 = DUP + DIT_preSess1.days 
#df$t1 <- df$DUP.months. + df$DIT_preSess1.days./30

#t1 = DIT (duration in treatment preSess1)
df$t1 <- df$DIT_preSess1.days./30

#time difference in "months" = weeks/4
df$t2<- df$t1 + as.numeric(difftime(df$MRI_2, df$MRI_1, units = 'weeks'))/4  
df$t3<- df$t1 + as.numeric(difftime(df$MRI_3, df$MRI_1, units = 'weeks'))/4
df$t4<- df$t1 + as.numeric(difftime(df$MRI_4, df$MRI_1, units = 'weeks'))/4

## Reorder by t1####
df <- df[order(df$t1), ]
row.names(df) <- NULL

df$plotID<-as.numeric(row.names(df))
specific_ticks <- c(12, 12*2, 12*4, 12*6, 12*7)
specific_labels <- c("","","","","")

## Plot ####
ggplot(df %>% 
         pivot_longer(cols=starts_with("t"),
                      names_to = c(".value", "sess"),
                      names_pattern = "(.+)(.)") %>%
         drop_na(t) %>% 
         mutate(sess = as.numeric(sess)), 
       aes(x = t)) +
  geom_line(aes(y=reorder(ID, plotID), group = ID ), size = .3) +
  geom_point(aes(y=reorder(ID, plotID), color = as.factor(Grupo),shape = factor(sess)), size=2) +
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


# Export the data frame to a CSV file
long_data <- pivot_longer(df, 
                          cols = starts_with("MRI") | starts_with("t"),
                          names_to = c(".value", "session"),
                          names_pattern = "(\\D+)(\\d+)")

write.csv(long_data, file = "subjects_finallist_DIT.csv", row.names = FALSE)



# PLOT 1b: DURATION OF ILLNESS####

## Exclude if not DUP ####
df <- df[!is.na(df$DUP.months.), ]

## Define time-points for plot ####

# t0 = 0 : inicio de enfermedad
# t1: time from illness start to first MRI scan
#t1 = DUP + DIT_preSess1.days 
df$t1 <- df$DUP.months. + df$DIT_preSess1.days./30

#time difference in "months" = weeks/4
df$t2<- df$t1 + as.numeric(difftime(df$MRI_2, df$MRI_1, units = 'weeks'))/4  
df$t3<- df$t1 + as.numeric(difftime(df$MRI_3, df$MRI_1, units = 'weeks'))/4
df$t4<- df$t1 + as.numeric(difftime(df$MRI_4, df$MRI_1, units = 'weeks'))/4

## Reorder by t1####
df <- df[order(df$t1), ]
row.names(df) <- NULL

df$plotID<-as.numeric(row.names(df))
specific_ticks <- c(12, 12*2, 12*4, 12*6, 12*7)
specific_labels <- c("","","","","")

## Plot ####
ggplot(df %>% 
         pivot_longer(cols=starts_with("t"),
                      names_to = c(".value", "sess"),
                      names_pattern = "(.+)(.)") %>%
         drop_na(t) %>% 
         mutate(sess = as.numeric(sess)), 
       aes(x = t)) +
  geom_line(aes(y=reorder(ID, plotID), group = ID ), size = .3) +
  geom_point(aes(y=reorder(ID, plotID), color = as.factor(Grupo),shape = factor(sess)), size=2) +
  labs(x = "Duration of illness (months)", y = "Subject ID", color = "Group", shape= "Session") +
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

# Export the data frame to a CSV file

long_data <- pivot_longer(df, 
                          cols = starts_with("MRI") | starts_with("t"),
                          names_to = c(".value", "session"),
                          names_pattern = "(\\D+)(\\d+)")

write.csv(df, file = "subjects_finallist_DUP.csv", row.names = FALSE)


