# Description ####
# Longitudinal plots to characterize our data

# Set environment ####
rm(list = ls()) # ctrl + L to clear console
setwd("/Volumes/PortableSSD/PostdocUC_projects/FEP_Longitudinal/")
library(tidyverse)


# PREPARE DATA ####
## Import data ####  
sample<-read.csv("patients_longitudinal.csv")
sample_HC<-read.csv("controls_longitudinal.csv")

#keep only columns and rows of interest
df <- sample[, c("ID", "Grupo","Fecha_nacimiento","Sexo", 
                          "MRI_1","MRI_2","MRI_3","MRI_4",
                          "DUP.months.","DIT_preSess1.days.",
                          "DIT_sess1toLastsess","Duration_of_illness.months.",
                          "Excluir")]
df <- df[1:117,] #ROWS up to subject 1517
rm(sample) #remove original full df

#same for controls
df_HC <- sample_HC[, c("ID", "Grupo","Fecha_nacimiento","Sexo", 
                 "MRI_1","MRI_2","MRI_3","MRI_4",
                 "Excluir")]
df_HC <- df_HC[1:88,] #ROWS up to subject C87

rm(sample_HC)


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

### CONTROLS####
#Correct class and format of dates (Birth and MRI dates)
df_HC$Fecha_nacimiento <- as.Date(df_HC$Fecha_nacimiento, format = "%d/%m/%Y")
df_HC$MRI_1 <- as.Date(df_HC$MRI_1, format = "%d/%m/%Y")
df_HC$MRI_2 <- as.Date(df_HC$MRI_2, format = "%d/%m/%Y")
df_HC$MRI_3 <- as.Date(df_HC$MRI_3, format = "%d/%m/%Y")
df_HC$MRI_4 <- as.Date(df_HC$MRI_4, format = "%d/%m/%Y")

#Factors: ID, Group (TRT, nonTRT), Sex
df_HC$ID <- as.factor(df_HC$ID)
df_HC$Grupo <- as.factor(df_HC$Grupo)
df_HC$Sexo <- as.factor(df_HC$Sexo)

sapply(df_HC, class)

## Exclude and remove NA ####
#PATIENTS
df <- df[!df$Excluir==1,]

df <- df[!is.na(df$Fecha_nacimiento), ]
df <- df[!is.na(df$Sexo), ]
df <- df[!is.na(df$DUP.months.), ]
df <- df[!is.na(df$DIT_preSess1.days.), ]

#CONTROLS
df_HC <- df_HC[!df_HC$Excluir==1,]

df_HC <- df_HC[!is.na(df_HC$Fecha_nacimiento), ]
df_HC <- df_HC[!is.na(df_HC$Sexo), ]


## Correct DATA ERRORS ####
df[df$ID=="1439","Fecha_nacimiento"] <- as.Date("11/07/1996",format = "%d/%m/%Y" )


#esta es la opción de incluirlos haciendo cero los NA
#df$DUP.months.[is.na(df$DUP.months.)] <- 0  



# PLOT 1: Duration of illness ####
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
specific_ticks <- c(12, 12*2, 12*4, 12*6, 12*8)
specific_labels <- c("","","","","")

## Plot####
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
  geom_text(aes(x=12*8, y=-20, label="8y"), size=3)




# PLOT 2: Age : ONLY PATIENTS####
## Define time-points for plot ####

# age0 = 0  : birth
# age1: age at sess1
df$age1 <- as.numeric(difftime(df$MRI_1, df$Fecha_nacimiento, units = 'weeks'))/52  

#time difference in "years" = weeks/52
df$age2 <-df$age1 + as.numeric(difftime(df$MRI_2, df$MRI_1, units = 'weeks'))/52
df$age3 <-df$age1 + as.numeric(difftime(df$MRI_3, df$MRI_1, units = 'weeks'))/52
df$age4 <-df$age1 + as.numeric(difftime(df$MRI_4, df$MRI_1, units = 'weeks'))/52

## Reorder by age1####
df <- df[order(df$age1), ]
row.names(df) <- NULL
df$plotID<-as.numeric(row.names(df))

## Plot####
ggplot(df %>% 
         pivot_longer(cols=starts_with("age"),
                      names_to = c(".value", "sess"),
                      names_pattern = "(.+)(.)") %>%
         drop_na(age) %>% 
         mutate(sess = as.numeric(sess)), 
       aes(x = age)) +
  geom_line(aes(y=reorder(ID, plotID), group = ID ), size = .3) +
  geom_point(aes(y=reorder(ID, plotID), color = as.factor(Grupo),shape = factor(sess)), size=2) +
  labs(x = "Age (y)", y = "Subject ID", color = "Group", shape= "Session") +
  theme(legend.position = "top", axis.text.y = element_text(size = 5)) +
  geom_density(aes(y = -..density..*100, group = Grupo, fill = Grupo, alpha=0.5)) +
  guides(fill = "none", alpha = "none") +
  geom_label(aes(x=33, y=-5, label="density"))




# PLOT 3: Age : ONLY CONTROLS####
## Define time-points for plot ####

# age0 = 0  : birth
# age1: age at sess1
df_HC$age1 <- as.numeric(difftime(df_HC$MRI_1, df_HC$Fecha_nacimiento, units = 'weeks'))/52  

#time difference in "years" = weeks/52
df_HC$age2 <-df_HC$age1 + as.numeric(difftime(df_HC$MRI_2, df_HC$MRI_1, units = 'weeks'))/52
df_HC$age3 <-df_HC$age1 + as.numeric(difftime(df_HC$MRI_3, df_HC$MRI_1, units = 'weeks'))/52
df_HC$age4 <-df_HC$age1 + as.numeric(difftime(df_HC$MRI_4, df_HC$MRI_1, units = 'weeks'))/52

## Reorder by age1####
df_HC <- df_HC[order(df_HC$age1), ]
row.names(df_HC) <- NULL
df_HC$plotID<-as.numeric(row.names(df_HC))

## Plot####
ggplot(df_HC %>% 
         pivot_longer(cols=starts_with("age"),
                      names_to = c(".value", "sess"),
                      names_pattern = "(.+)(.)") %>%
         drop_na(age) %>% 
         mutate(sess = as.numeric(sess)), 
       aes(x = age)) +
  geom_line(aes(y=reorder(ID, plotID), group = ID ), size = .3) +
  geom_point(aes(y=reorder(ID, plotID), color = as.factor(Grupo),shape = factor(sess)), size=2) +
  labs(x = "Age (y)", y = "Subject ID", color = "Group", shape= "Session") +
  theme(legend.position = "top", axis.text.y = element_text(size = 5)) +
  geom_density(aes(y = -..density..*100, group = Grupo, fill = Grupo, alpha=0.5)) +
  guides(fill = "none", alpha = "none") +
  geom_label(aes(x=33, y=-5, label="density"))






# PLOT 4: Age : FULL DATA SAMPLE####

df_full <- rbind(df[,1:8], df_HC[1:8])

## Define time-points for plot ####

# age0 = 0  : birth
# age1: age at sess1
df_full$age1 <- as.numeric(difftime(df_full$MRI_1, df_full$Fecha_nacimiento, units = 'weeks'))/52  

#time difference in "years" = weeks/52
df_full$age2 <-df_full$age1 + as.numeric(difftime(df_full$MRI_2, df_full$MRI_1, units = 'weeks'))/52
df_full$age3 <-df_full$age1 + as.numeric(difftime(df_full$MRI_3, df_full$MRI_1, units = 'weeks'))/52
df_full$age4 <-df_full$age1 + as.numeric(difftime(df_full$MRI_4, df_full$MRI_1, units = 'weeks'))/52

## Reorder by age1####
df_full <- df_full[order(df_full$age1), ]
row.names(df_full) <- NULL
df_full$plotID<-as.numeric(row.names(df_full))

## Plot####
ggplot(df_full %>% 
         pivot_longer(cols=starts_with("age"),
                      names_to = c(".value", "sess"),
                      names_pattern = "(.+)(.)") %>%
         drop_na(age) %>% 
         mutate(sess = as.numeric(sess)), 
       aes(x = age)) +
  geom_line(aes(y=reorder(ID, plotID), group = ID ), size = .3) +
  geom_point(aes(y=reorder(ID, plotID), color = as.factor(Grupo),shape = factor(sess)), size=2) +
  labs(x = "Age (y)", y = "Subject ID", color = "Group", shape= "Session") +
  theme(legend.position = "top", axis.text.y = element_text(size = 5)) +
  geom_density(aes(y = -..density..*100, group = Grupo, fill = Grupo, alpha=0.5)) +
  guides(fill = "none", alpha = "none") +
  geom_label(aes(x=33, y=-5, label="density"))

