###Zapata Luis, Caravagna Giulio 2019
# ###Simulation of Population growth, dNdS and mutation accumulation based on probabilities 
library(ggpubr)
library(tidyverse)
library(ggmuller)
library(crayon)
library(pio)
library(easypar)
library(parallel)

rm(list = ls())
source('simulator.R')
source("functions_ZCsimulator.R")

set.seed(12345)
######### Number of repetitions #########
TASKS = 50

######## Define our set of parameters for the SIESTA model ###########
params = NULL

######## Global parameters #########
params['GENERATIONS'] = 100
params['Kcapacity'] = 2000

#######Set probability of cell survival
params['ps'] <- 0.5
params['pdiff'] <- (1 - params['ps'])

######Set mutation occurrence parameters
params['u'] = 1 /  (10 ^ 6) #muts/ p divisionm / p bp, ##mutation occurrence/ polymerase error
params['L'] = 50 * (10 ^ 6) ###Length of the genome/exome size of exome syn plus nonsyn mutations
params['C'] = 1 /  (10 ^ 2) #correction efficiency, it allows only 1 out of X number of mutations to be preserved at each cell division

####Set value for the probability of the immune system to kill an immungenic population
params['pattack'] = 0

#######Set INITIAL conditions everytime yuo restart have to uopdate all probabilities
params['pnad'] <- 0.00 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.00 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.000 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations



####################################################################################################################
################################ New simulation ####################################################################


params['u'] = 1 / (10 ^ 6)
#######Set mutation probabilities
params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.00 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.000 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN4 = lapply(1:TASKS,function(w)  list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN4,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_01.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################


params['u'] = 1 / (10 ^ 6)
#######Set mutation probabilities
params['pnad'] <- 0.005 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.00 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.000 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN4 = lapply(1:TASKS,function(w)  list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN4,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_005.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################


params['u'] = 1 / (10 ^ 6)
#######Set mutation probabilities
params['pnad'] <- 0.001 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.00 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.000 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN4 = lapply(1:TASKS,function(w)  list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN4,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_001.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver  plus 5% immunogenic + pattack 0%

#######Set mutation probabilities
params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.05 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.000 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations
params['pattack'] <- 0

PARAMS_RUN6 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN6,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_00_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver plus 5% immunogenic + pattack 5%

params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.05 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.00 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations
params['pattack'] <- 0.05

PARAMS_RUN7 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN7,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_05_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver plus deleterious plus 5% immunogenic + pattack 25% 

params['pattack'] <- 0.25

PARAMS_RUN8 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN8,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_25_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver plus deleterious plus 5% immunogenic + pattack 50% 

params['pattack'] <- 0.5

PARAMS_RUN9 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN9,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_50_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver plus deleterious plus 5% immunogenic + pattack 50% + 5% escape

params['pattack'] <- 0.75

PARAMS_RUN80 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN80,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_75_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver plus deleterious plus 5% immunogenic + pattack 50% + 5% escape

params['pattack'] <- 0.90

PARAMS_RUN81 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN81,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_90_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################
###Driver plus deleterious plus 5% immunogenic + pattack 50% + 5% escape

params['pattack'] <- 0.99

PARAMS_RUN82 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN82,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_99_driver_01_IM_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################

params['pattack'] <- 0.99

#######Set INITIAL conditions everytime yuo restart have to uopdate all probabilities
params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.05 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.01 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN100 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN100,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_99_driver_01_IM_05_E_01.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################

params['pattack'] <- 0.99

#######Set INITIAL conditions everytime yuo restart have to uopdate all probabilities
params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.05 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.05 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN101 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN101,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_99_driver_01_IM_05_E_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################

params['pattack'] <- 0.99

#######Set INITIAL conditions everytime yuo restart have to uopdate all probabilities
params['pnad'] <- 0.001 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.05 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.05 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN102 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN102,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_99_driver_001_IM_05_E_05.RData")

####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################

params['pattack'] <- 0.99

#######Set INITIAL conditions everytime yuo restart have to uopdate all probabilities
params['pnad'] <- 0.001 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.02 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.01 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN103 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN103,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_99_driver_001_IM_02_E_01.RData")


####################################################################################################################
################################ New simulation ####################################################################
####################################################################################################################

params['pattack'] <- 0.99

#######Set INITIAL conditions everytime yuo restart have to uopdate all probabilities
params['pnad'] <- 0.001 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.02 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.05 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

PARAMS_RUN104 = lapply(1:TASKS,function(w)  list(params))
RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN104,
  silent = TRUE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_IA_99_driver_001_IM_02_E_05.RData")