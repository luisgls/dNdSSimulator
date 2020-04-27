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
source("~/tools/dNdSSimulator/simulator_modelA.R")
source("~/tools/dNdSSimulator/functions_ZCsimulator.R")

#set.seed(12345)
dummy<-as.integer(runif(1,min=1,max=100000))
set.seed(dummy)

######### Number of repetitions #########
TASKS = 4

######## Define our set of parameters for the SIESTA model ###########
params = NULL
####Define initial population size
params['n0'] <- 5

######## Global parameters #########
params['GENERATIONS'] = 50
params['Kcapacity'] = 1000

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
  parallel = FALSE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)

RUNS.df = plyr::ldply(RUNS, as.data.frame)
write.csv(RUNS.df, file = paste(dummy,"simulation_TABLE.csv",sep="_"), row.names = F, quote = F)

#save(RUNS, file=paste(dummy,"simulation.RData",sep="."))
all_genotypes<-lapply(RUNS,get_sequencing)
SEQ.df = plyr::ldply(all_genotypes, as.data.frame)
write.csv(SEQ.df, file = paste(dummy,"simulation_SEQ.csv",sep="_"), row.names = F, quote = F)

#save(all_genotypes, file=paste(dummy,"simulation_SEQ.RData",sep="."))

# 
# PARAMS_DT<-Reduce(bind_rows, lapply(PARAMS_RUN1, data.frame))
# lapply(PARAMS_RUN1, data.frame)
# 
# # PAR = Reduce(bind_rows, lapply(PAR, data.frame))
# # TASKS = data.frame(
# #   pnad = c(.3, .23),
# #   pnac = c(.4, .2),
# #   pnae = c(.8, .23)
# # )
# out_logs = 'logs'
# task = 'mytask'
# # CREAT RUNNING SCRIPT
# dir.create(out_logs)
# cluster_config = easypar::default_BSUB_config()
# cluster_config$`-J` = "LUIS"
# cluster_config$`-P` = ""
# cluster_config$`-q` = "normal"
# cluster_config$`-o` = paste0('./', out_logs, '/out.%J.%I')
# cluster_config$`-e` = paste0('./', out_logs, '/err.%J.%I')
# # Run
# easypar::run_lsf(
#   FUN = simulator,
#   PARAMS = PARAMS_DT,
#   extra_commands = 'R_LIBS_USER="-----"',
#   BSUB_config = cluster_config,
#   modules = 'R/3.6.0',
#   run = FALSE,
#   Submission_script = paste0(task, ".sh"),a
#   input_file = paste0(task, "_input_jobarray.csv"),
#   R_script = paste0(task, "_rscript.R")
# )
# 
# RUNS = run(
#   FUN = simulator,
#   PARAMS = PARAMS_RUN1,
#   silent = FALSE,
#   parallel = TRUE,
#   packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# Message caravagn
# 
