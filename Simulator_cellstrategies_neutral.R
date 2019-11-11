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
TASKS = 100

######## Define our set of parameters for the SIESTA model ###########
params = NULL

###Initial population size
params['n0'] = 5 ## 2^5 = 32 cells

#######Set mutation probabilities
params['pnad'] <- 0.00 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.00 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.000 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

#######Set probability of cell survival
params['ps'] <- 0.5
params['pdiff'] <- (1 - params['ps'] )

######Set mutation occurrence parameters
params['u'] = 1 /  (10 ^ 6) #muts/ p divisionm / p bp, ##mutation occurrence/ polymerase error
params['L'] = 50 * (10 ^ 6) ###Length of the genome/exome size of exome syn plus nonsyn mutations
params['C'] = 1 /  (10 ^ 2) #correction efficiency, it allows only 1 out of X number of mutations to be preserved at each cell division


####Set value for te probability of te immune system to kill an immungenic population
params['pattack'] = 0
params['GENERATIONS'] = 100
params['Kcapacity'] = 2000

PARAMS_RUN1 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

# run(
#   FUN = plotter,
#   PARAMS = lapply(1:TASKS, function(w) list(j = w, M = RUNS[[w]])),
#   silent = FALSE,
#   parallel = FALSE,
#   packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
# )

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations/neutral_MSS.RData")

###Mid Mutation rate
params['u'] = 10 / (10 ^ 6)
PARAMS_RUN2 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN2,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations/neutral_MSI.RData")

###High Mutation rate
params['u'] = 100 / (10 ^ 6)
PARAMS_RUN3 = lapply(1:TASKS,function(w)  list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN3,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations/neutral_POLE.RData")


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
TASKS = 100

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


#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral.RData")

#NEUTRAL MSI, MID MUTATRION RATE
params['u'] = 10 / (10 ^ 6)

PARAMS_RUN2 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN2,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSI_neutral.RData")

#NEUTRAL POLE, HIGH MUTATRION RATE
params['u'] = 100 / (10 ^ 6)

PARAMS_RUN3 = lapply(1:TASKS,function(w)  list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN3,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/POLE_neutral.RData")

###################################################################################################
###################Run for different initial population sizes######################################
###################################################################################################

params['u'] = 1 /  (10 ^ 6) #muts/ p divisionm / p bp, ##mutation occurrence/ polymerase error

params['n0'] = 3 
#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1001 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1001,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral_8.RData")


params['n0'] = 4
#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1002 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1002,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral_16.RData")

params['n0'] = 5
#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1003 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1003,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral_32.RData")

params['n0'] = 6
#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1004 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1004,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral_64.RData")

params['n0'] = 7
#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1005 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1005,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral_128.RData")

params['n0'] = 8
#NEUTRAL MSS, LOW MUTATRION RATE
PARAMS_RUN1006 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1006,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

save(RUNS, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_neutral_256.RData")