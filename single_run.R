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

set.seed(NULL)
######### Number of repetitions #########
TASKS = 1

######## Define our set of parameters for the SIESTA model ###########
params = NULL

####Define initial population size
params['n0'] <- 5

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

#######Set mutation probabilities
params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.10 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.01 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

####Set value for the probability of the immune system to kill an immungenic population
params['pattack'] = 0.5


PARAMS_RUN1 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

run(
  FUN = plotter,
  PARAMS = lapply(1:TASKS, function(w) list(j = w, M = RUNS[[w]])),
  silent = FALSE,
  parallel = FALSE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
