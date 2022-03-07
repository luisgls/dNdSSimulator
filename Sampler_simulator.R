library(ggpubr)
library(tidyverse)
library(crayon)
library(pio)
library(easypar)
library(parallel)
source('simulator_modelA.R')
source("functions_ZCsimulator.R")
source("simulator_plotting.R")

######## Grid of values

rangeify = function(x, y, n = 10) seq(x, y, (y-x)/n)

# mutation rate, etc.
u_range = c(1e-04, 1e-06, 1e-08)
pnad_range = rangeify(0.001, 0.05)
pnai_range = rangeify(0.01, 0.2)
pnae_range = rangeify(0.001, 0.05)
pattack_range = rangeify(0.01, 1)
repeats = 1:10

grid = expand_grid(u_range, pnad_range, pnai_range, pnae_range, pattack_range, r = repeats)

saveRDS(grid, 'grid.RDS')

runner = function(i)
{
  grid = readRDS('grid.RDS')
  cli::cli_h3("Test {.field {i}} / {.value {nrow(grid)}}")
  TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

  ######## Define a Null set of parameters for our evolutionary model ###########
  params = NULL

  ####Define initial population size (number of wild-type cells)
  params['n0'] <- 5

  ######## Global parameter to stop simulation after N generations #########
  params['GENERATIONS'] = 75

  ######## Global parameter to stop simulation after K number of alive cells #########
  params['Kcapacity'] = 1000

  ####### Set probability of cell survival #########
  params['ps'] <- 0.5
  params['pdiff'] <- (1 - params['ps'])

  ###### Set mutation rate parameter, the final mutation rate value used is equal to a combination of the following terms: u*C
  params['u'] = grid$u_range[i] ### muts/ division / bp, mutation occurrence/ polymerase error
  params['L'] = 50 * (10 ^ 6) ### Length of the genome/exome size of exome syn plus nonsyn mutations
  params['C'] = 1 /  (10 ^ 2) ### damage repair correction efficiency, it allows only 1 out of 100 number of errors to be preserved at each cell division

  #######Set mutation probabilities
  params['pnad'] <- grid$pnad_range[i] #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate, etc)
  params['pnai'] <- grid$pnai_range[i] #Nonsyn - deleterious mutations (immunogenic mutations)
  params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
  params['pnae'] <- grid$pnae_range[i] #Nonsyn - escape immune system mutations
  params['pnsi'] <- params['pnai']/3    # Synonymous mutations in immunopeptidome
  params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
  params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
  params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

  s_one = params['pnad'] + params['pnai'] +  params['pnae'] +  params['pnak'] +   params['pnap'] +  params['pns']
  if(s_one != 1) {
    params['pnad'] = params['pnad']/s_one
    params['pnai'] = params['pnai']/s_one
    params['pnae'] = params['pnae']/s_one
    params['pnak'] = params['pnak']/s_one
    params['pnap'] = params['pnap']/s_one
    params['pns'] = params['pns']/s_one
  }

  ####Set value for the probability of the immune system to kill an immungenic population (PIs in the manuscript)
  params['pattack'] = grid$pattack_range[i]

  M<-simulator(params, id = i)

  M = list(input = params, output = M)

  saveRDS(M, paste0(i, '.RDS'))

  TIME = difftime(as.POSIXct(Sys.time(), format = "%H:%M:%S"),
                  TIME, units = "mins")
  TIME = round(TIME, 2)
  cli::cli_h3("Test {.field {i}} / {.value {nrow(grid)}} COMPLETED - {.value {TIME}} (mins) ")

}

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)



