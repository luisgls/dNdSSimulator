
# dNdSSimulator
dNdS Simulator uses R to simulate a stochastic branching process with a set of genotypes and phenotypes per cell. The main advantage of the model is that simulates immunoediting by introducing immunogenic and escape mutations as well as driver and passenger mutations.


## Installation dN/dS simulation package

We must first clone our repository, load libraries and functions for the simulation:

```{bash clone,eval=F}
git clone https://github.com/luisgls/dNdSSimulator.git
```

```{r dependencies, message=F}
library(ggpubr)
library(tidyverse)
library(crayon)
library(pio)
library(easypar)
library(parallel)
source('simulator_modelA.R')
source("functions_ZCsimulator.R")
source("simulator_plotting.R")
```
*Note*:
Libraries pio and easypar can be obtained from <https://caravagn.github.io/>

*Note*:
If you want to run the second model where the immune system attacks a population of immunogenic cells, you must source simulator_modelB.R instead of simulator_modelA.R.

## Running one simulation

As an example, we will run the simulator using the following paramters:

```{r parameters}

######### Number of simulations that we will perform #########
TASKS = 1

######## Define a Null set of parameters for our evolutionary model ###########
params = NULL

####Define initial population size (number of wild-type cells)
params['n0'] <- 5

######## Global parameter to stop simulation after N generations #########
params['GENERATIONS'] = 50

######## Global parameter to stop simulation after K number of alive cells #########
params['Kcapacity'] = 500

####### Set probability of cell survival #########
params['ps'] <- 0.5
params['pdiff'] <- (1 - params['ps'])

###### Set mutation rate parameter, the final mutation rate value used is equal to a combination of the following terms: u*C
params['u'] = 1 /  (10 ^ 6) ### muts/ division / bp, mutation occurrence/ polymerase error
params['L'] = 50 * (10 ^ 6) ### Length of the genome/exome size of exome syn plus nonsyn mutations
params['C'] = 1 /  (10 ^ 2) ### damage repair correction efficiency, it allows only 1 out of 100 number of errors to be preserved at each cell division

#######Set mutation probabilities
params['pnad'] <- 0.01 #Nonsyn - driver mutations (increase proliferation rate, oncogenes/ reduce correction efficiency, mut_modifiers, decrease death rate, etc)
params['pnai'] <- 0.10 #Nonsyn - deleterious mutations (immunogenic mutations)
params['pnak'] <- 0.00 #Nonsyn - deleterious mutations (death)
params['pnae'] <- 0.01 #Nonsyn - escape immune system mutations
params['pnsi'] <- params['pnai']/3    # Synonymous mutations in immunopeptidome
params['pnsd'] <- params['pnad']/3    # Synonymous mutations in driver
params['pnap'] <- 0.75 - ( params['pnad'] +  params['pnai'] +  params['pnae'] +  params['pnak']) #Nonsyn - passenger mutations
params['pns']  <- 0.25 - ( params['pnsi'] +  params['pnsd']) #Synonymous mutations

####Set value for the probability of the immune system to kill an immungenic population (PIs in the manuscript)
params['pattack'] = 0.5
```

Now that we have defined the parameters for the simulation, we will run one simulation using the command simulator:

```{r simul1, cache=T}
M<-simulator(params)
```

The output contains the following information for each single cell:

```{r eval=F}
id  #cell id
pnad, pnap, pnai, pnae, pnak, pns, pnsi, pnsd #vector of probabilities to acquire a genotype
gnad, gnap, gnai, gnae, gnak, gns, gnsi, gnsd #vector of genotype
strategy # Death (never replicated), Survive (Replicate), immunogenic or escape cell
parent # parental cell ID
alive # status at the end of simulation
t_0 # current time
t_e # effective time
dnds # global dN/dS of the cell
dnds_i # immune dN/dS of the cell 
dnds_d # driver dN/dS of the cel
psurv # Probability of effective cell division (delta)
pimmune # Probability of immune recognition
clone_id # Clone id assigned by driver status
clone_I_id # Clone id assigned by immunogenic status
parent_clone_id # Parental clone ID (driver status)
```
## Running multiple simulations

We have implemented the following code to run multiple simulations in parallel

```{r simul2, cache=T, eval=F}
######### Number of simulations that we will perform #########
TASKS = 1000

PARAMS_RUN1 = lapply(1:TASKS,function(w) list(params))

RUNS = run(
  FUN = simulator,
  PARAMS = PARAMS_RUN1,
  silent = FALSE,
  parallel = TRUE,
  packages = c("ggpubr","tidyverse","ggmuller","crayon","pio")
)
```
This will create a list of objects, each object corresponds to one simulation.

## Visualization

We  can also plot the results using the following function (after loading some libraries):

```{r visualize,message=F,warning=F}
###Libraries for visualizing
library(reshape2)
library(tidygraph)
library(ggraph)
library(ggstatsplot)
library(ggmuller)

###Code to visualize results
  # Visualization
  if (count(M %>% dplyr::select(clone_id) %>% unique()) > 1) {
    p1 <-
      mullerplot(
        M,
        t_from = 1,
        #  t_to = 50,
        legend = TRUE,
        ps_0 = ps,
        cutoff = .05,
        add_clone_tree = FALSE
      )
  }

  p1
```
which produces the following figure:

https://github.com/luisgls/dNdSSimulator/blob/master/dNdSTutorial_files/figure-html/visualize-1.png

https://htmlpreview.github.io/?https://github.com/luisgls/dNdSSimulator/blob/master/dNdSTutorial.html

## Citation

If you want to use the simulator please cite the following article: 

<https://www.biorxiv.org/content/10.1101/2020.07.21.215038v2>
