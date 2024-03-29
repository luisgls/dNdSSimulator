###Zapata Luis, Caravagna Giulio 2019
# ###Simulation of Population growth, dNdS and mutation accumulation based on probabilities 

##Main simulator function to repeat multiple times given a set of paramteres
simulator = function(list_params)
  {
  library(ggpubr)
  library(tidyverse)
  library(ggmuller)
  library(crayon)
  library(pio)
  library(easypar)
  library(parallel)
  
  source("functions_ZCsimulator.R")
  ###Initital number of cell
  n0 = list_params['n0']
  
  ###initial genotypes
  pnad = list_params['pnad']
  pnai = list_params['pnai']
  pnak = list_params['pnak']
  pnae = list_params['pnae'] 
  pnsi = list_params['pnsi'] 
  pnsd = list_params['pnsd'] 
  pnap = list_params['pnap'] 
  pns = list_params['pns']  
  
  #######Set probability of cell survival
  ps = list_params['ps'] 
  pdiff = list_params['pdiff']
  
  ######Set mutation occurrence parameters
  u = list_params['u']
  L = list_params['L']  
  C = list_params['C']
  
  
  ####Set value for the probability of the immune system to kill an immunogenic population
  pattack = list_params['pattack']
  GENERATIONS = list_params['GENERATIONS']
  Kcapacity = list_params['Kcapacity']
  
  N = u * L * C ##The effective average number of coding mutations per each cell division

  #Vector of types of mutations
  a0 <- c("gnad", "gnap", "gnai", "gnae", "gnak", "gns", "gnsi", "gnsd")
  p0 <- c(pnad , pnap , pnai , pnae , pnak , pns , pnsi, pnsd)
  
  #Initialize empty data frame
  M <- initial_state(as.list(p0), ps, n0)
  
  ####Program
  #pb = txtProgressBar(min = 0, max = GENERATIONS, style = 3)
  
  #####Store bottlenecks
  bottleneck_time = NA
  
  #set.seed(1234567)
  set.seed(NULL)
  
  # create new daughter cell
  for (i in 1:GENERATIONS) {
    
    ##remove from the memory death cells
    #M <- M %>% filter(strategy != "Death") 

    ##Get alive cells to d+o the looping
    A = get_alive(M)
    
    ##Carrying capacity of the population
    if (length(A) > Kcapacity) {
      break
      #M <- M %>% mutate(psurv = ifelse(psurv > 0.5, 0.5, psurv))
    }
    
    ###For each cell alive in the previous generation creates two daugther cells and kill the parent if they are alive
    for (cell in A) {
      ##Kill the cell if it has any immunogenic mutation
      #M = set_immunedead1(M, cell, i, pattack)
      
      ##Kill the cell if it is immunogenic (immune) and the immune system is active (Model 2. Kill immunogenic cells based on probabilities)
      M = set_immunedead2(M, cell, i, pattack)
    
      #Start cell division and copying genotype from parental cell
      if (get_divide(M, cell) & is_alive(M, cell)) {
        ##Cell divides two times
        M = simulate_daughter_cell(M, cell, i, a0, N)
        M = simulate_daughter_cell(M, cell, i, a0, N)
        
        ##Kill paternal cell
        M = set_dead_cell(M, cell, i)
        
      }
      else{
        ###Cells do not go to division
        M = kill_cell(M, cell, i)
      }
    }
    
  }
  return(M)
}
