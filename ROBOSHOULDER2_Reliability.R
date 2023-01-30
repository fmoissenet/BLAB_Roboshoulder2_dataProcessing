# ROBOSHOULDER2 Project
# Compute reliability and interpretability parameters (ICC, SEM, MDC)
# with a Fit Linear Mixed-Effects Models

# Clear workspace
rm(list = ls())

# Load libraries 
library(rstudioapi)
library(tidyr)
library(magrittr)
library(lme4)
library(ICC.Sample.Size)

# Select working directories
folder_inputs  <- "C:/Users/moissene/OneDrive - unige.ch/Article ROBOSHOULDER2/Données/Stats/input"
folder_outputs <- "C:/Users/moissene/OneDrive - unige.ch/Article ROBOSHOULDER2/Données/Stats/output"
setwd(folder_inputs)

# Set dataset 
Parameters <- list.files(folder_inputs)
tmp        <- dim.data.frame(Parameters)
  
# Initialisation
Unit       <- 0
ICC_intra  <- 0
ICC_inter  <- 0
SEM_intra  <- 0
SEM_inter  <- 0
SEMp_intra <- 0
SEMp_inter <- 0
MDC_intra  <- 0
MDC_inter  <- 0
MDCp_intra <- 0
MDCp_inter <- 0
  
# START - Loop on parameters
for (i in 1:tmp[2]){
  # Load data
  file       <- Parameters[i]
  Input_Data <- read.csv(paste(folder_inputs,file,sep="/"),header=TRUE,sep=",",dec=".")
  Unit[i]    <- Input_Data$UNIT[1]  
  dim(Input_Data)
  Input_Data <- na.omit(Input_Data)
    
  # Prepare data
  Input_Data$CARTILAGE      <- factor(Input_Data$CARTILAGE, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9")) 
  Input_Data$OPERATOR       <- factor(Input_Data$OPERATOR, levels=c("1", "2", "3"))
  Input_Data$DIGITALISATION <- factor(Input_Data$DIGITALISATION, levels=c("1", "2", "3"))   
     
  # Compute variance Components with LMER Model (Fit Linear Mixed-Effects Models)
  mod       <- lmer(VALUE~1 + (1|CARTILAGE) + (1|CARTILAGE:OPERATOR) + (1|CARTILAGE:DIGITALISATION), data = Input_Data)
  variance  <- as.data.frame(VarCorr(mod))
  v_car_dig <- variance$vcov[1] # Variance associated with Digitalisation
  v_car_ope <- variance$vcov[2] # Variance associated with Operator
  v_car     <- variance$vcov[3] # Variance associated with Cartilage
  v_res     <- variance$vcov[4] # Residual Variance
      
  # Total variance
  v_tot <- v_car_dig + v_car_ope + v_car + v_res
      
  # Variance outcome
  v_intra <- v_car_dig + v_res
  v_inter <- v_car_ope + v_res
      
  # ICC
  ICC_intra[i] <- (v_tot - v_intra) / v_tot
  ICC_inter[i] <- (v_tot - v_inter) / v_tot
      
  # SEM
  SEM_intra[i] <- sqrt(v_tot * (1 - ICC_intra[i]))
  SEM_inter[i] <- sqrt(v_tot * (1 - ICC_inter[i]))
  
  # SEM%
  SEMp_intra[i] <- SEM_intra[i]/mean(Input_Data$VALUE[])*100
  SEMp_inter[i] <- SEM_inter[i]/mean(Input_Data$VALUE[])*100
      
  # MDC95
  MDC_intra[i] <- 1.96 * sqrt(2) * SEM_intra[i]
  MDC_inter[i] <- 1.96 * sqrt(2) * SEM_inter[i]
  
  # MDC%
  MDCp_intra[i] <- MDC_intra[i]/mean(Input_Data$VALUE[])*100
  MDCp_inter[i] <- MDC_inter[i]/mean(Input_Data$VALUE[])*100

  # Clean data
  rm("Input_Data", "mod", "variance", "v_car_dig", "v_car_ope", 
     "v_car", "v_res", "v_tot", "v_intra", "v_inter")
}
# END - Loop on parameters
  
# Export Data
Reliability      <- cbind(Parameters, Unit, ICC_intra,  ICC_inter, SEM_intra,  SEM_inter, SEMp_intra,  SEMp_inter)
Interpretability <- cbind(Parameters, Unit, MDC_intra,  MDC_inter, MDCp_intra,  MDCp_inter)
write.csv(Reliability,paste(folder_outputs, "Reliability.csv", sep="/"))
write.csv(Interpretability,paste(folder_outputs, "Interpretability.csv", sep="/"))