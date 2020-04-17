library(tidyverse)
library(ctree)
library(ggsci)
library(ggrepel)
library(tidyverse)
library("scales")
library(ggpubr)
library(RColorBrewer)
rm(list = ls())
source("functions_ZCsimulator.R")
#source("simulator.R")
source("simulator_plotting.R")


###############################    Neutral


###############################    Driver Only
# load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_01.RData")
# #Perform sequencing (This step is very time consumig) and save to Rdata file
# all_genotypes<-lapply(RUNS,get_sequencing)
# save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_01_SEQ.RData")
# 
# 
# load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_005.RData")
# all_genotypes<-lapply(RUNS,get_sequencing)
# save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_005_SEQ.RData")
# 
# 
# load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_001.RData")
# all_genotypes<-lapply(RUNS,get_sequencing)
# save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations3/MSS_D_001_SEQ.RData")


###############################    Driver + Immune Model A
load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_00_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_00_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_25_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_25_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_50_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_50_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_75_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_75_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_100_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_100_driver_01_IM_05_SEQ.RData")

###############################    Driver + Immune Model B
load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_00_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_00_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_25_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_25_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_50_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_50_driver_01_IM_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_75_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_75_driver_01_IM_05_SEQ.RData")


load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_100_driver_01_IM_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_100_driver_01_IM_05_SEQ.RData")


###############################    Driver + Immune A + Escape 
load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_99_driver_01_IM_05_E_01.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_99_driver_01_IM_05_E_01_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_99_driver_01_IM_05_E_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_99_driver_01_IM_05_E_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MSS_IA_99_driver_01_IM_05_E_001.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelA/MMSS_IA_99_driver_01_IM_05_E_001_SEQ.RData")



###############################    Driver + Immune B + Escape 
load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_99_driver_01_IM_05_E_01.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_99_driver_01_IM_05_E_01_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_99_driver_01_IM_05_E_05.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_99_driver_01_IM_05_E_05_SEQ.RData")

load(file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MSS_IA_99_driver_01_IM_05_E_001.RData")
all_genotypes<-lapply(RUNS,get_sequencing)
save(all_genotypes, file="~/Dropbox (Personal)/Postdoc/Projects/simul_dnds/simulations_modelB/MMSS_IA_99_driver_01_IM_05_E_001_SEQ.RData")

