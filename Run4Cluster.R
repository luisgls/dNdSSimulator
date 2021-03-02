###Zapata Luis, Caravagna Giulio 2019
# ###Simulation of Population growth, dNdS and mutation accumulation based on probabilities 

.libPaths( c("/home/haem/lortiz/R/3.5.0/", .libPaths()))

library(ggpubr)
library(tidyverse)
library(ggmuller)
library(crayon)
library(pio)
library(easypar)
library(parallel)

rm(list = ls())
source("~/tools/dNdSSimulator/simulator_modelA.R")
source("~/tools/dNdSSimulator/functions_ZCsimulator.R")

#set.seed(12345)
dummy<-as.integer(runif(1,min=1,max=1000000))
set.seed(dummy)

######### Number of repetitions #########
#TASKS = 4

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
params['pnad'] <- 0.00 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate)
params['pnai'] <- 0.00 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.00 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in imunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

####Set value for the probability of the immune system to kill an immungenic population
params['pattack'] = 0.0

RUN<-simulator(params)
write.csv(RUN, file = paste(dummy,"simulation_TABLE.csv",sep="_"), row.names = F, quote = F)
all_genotypes<-get_sequencing(RUN)
write.csv(all_genotypes, file = paste(dummy,"simulation_SEQ.csv",sep="_"), row.names = F, quote = F)
