##functions for Zapata Caravagna simulator

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Create a new simualtion
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

##Create an empty data frame
initial_state = function(p0, ps, generations)
{
  numCells = 2^generations
  
  M = Reduce(
    bind_rows,
    lapply(1:numCells, function(w) empty_dataframe(p0, ps))
  )
  
  M$id = 1:numCells
  M
}

empty_dataframe <- function(p0, ps)
{
  # State of the system
  M <- data.frame(
    id = NA,
    pnad = NA,
    pnap = NA,
    pnai = NA,
    pnae = NA,
    pnak = NA,
    pns = NA,
    pnsi = NA,
    pnsd = NA,
    gnad = 0,
    gnap = 0,
    gnai = 0,
    gnae = 0,
    gnak = 0,
    gns = 0,
    gnsi = 0,
    gnsd = 0,
    strategy = NA,
    parent = 0,
    alive = TRUE,
    t_0 = 0,
    t_e = 0,
    dnds = NA,
    dnds_i = NA,
    dnds_d = NA,
    psurv = NA,
    pimmune = NA,
    clone_id = NA,
    clone_I_id = NA,
    parent_clone_id = NA,
    stringsAsFactors = FALSE
  ) %>% as_tibble()
  
  # Initialization
  M$id = 1
  M[, c('pnad', 'pnap', 'pnai', 'pnae', 'pnak', 'pns', 'pnsi', 'pnsd')] = p0
  M[, c('gnad', 'gnap', 'gnai', 'gnae', 'gnak', 'gns', 'gnsi', 'gnsd')] = 0
  M$strategy = "Survive"
  #M$gnai = 3/32
  #M$gnap = 3/32
  #M$gnsi = 1/32
  M$parent = 0
  M$alive = T
  M$t_0 = 0
  M$dnds = 1
  M$dnds_i = 1
  M$dnds_d = 1
  M$psurv = ps
  M$pimmune = 0
  M$clone_id = 0
  M$clone_I_id = 0
  M$parent_clone_id = 0
  
  M
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Getter functions for the state of the model
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

##Get probability vector
get_prob = function(M, id) {
  M %>% filter(id == !!id) %>%
    select(starts_with('pn')) %>%
    unlist()
}

##Get genotype
get_geno = function(M, id) {
  M %>% filter(id == !!id) %>%
    select(starts_with('gn')) %>%
    unlist()
}

##Get alive cells 
get_alive = function(M) {
  M %>% filter(alive) %>% pull(id)
}

##Get new id when we create a cell
get_new_id = function(M) {
  max(M$id) + 1
}

#Get if a cell divides or not 
get_divide = function(M, parent_id) {
  M %>% filter(id == !!parent_id) %>% pull(psurv) >= runif(1)
}

is_alive = function(M, parent_id) {
  M %>% filter(id == !!parent_id) %>% pull(alive)
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Setter functions for the state of the model
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Update a genotype and applies changes to the clonal
# structure of a tumour as well (changes to the probabilities)
playStrategy <- function(M, id, new_geno, parent_id) 
{
  outcome = chooseStrategy(M, id)
  old_parent_clone_id = M %>% filter(id == !!parent_id) %>% pull(clone_id)
  newclone_id = max(M$clone_id) + 1
  newclone_I_id = max(M$clone_I_id) + 1
  
  ######################################################################################################################################################
  ############################################# PHENOTYPE 1: Probability of survival ###################################################################
  #Function to modify the probability of survival of the clone
  ##Delta value for increasing or decreasing survival probability based on the number of driver events
  #delta_decrease = 0
  #delta_increase = 0
  
  delta_decrease = rgamma(1, 0.5, rate = 10)
  x <- M %>% filter(id == !!id) %>% pull(gnad)
  #Old
  #delta_increase = ( 1/ (1+exp(1)^(-0.7*x) ) )/2 
  #New(21-5-2019)
  delta_increase = gompertz(0.5,5,1,x)
  
  change_psurv = function(psurv) 
  {
    # deleterious reduce the probability of survival
    if(new_geno['gnak'] >= 1) 
    {
      newval = psurv - delta_decrease
      if(newval <= 0.001) newval = 0.001
      return(newval)
    }    
    # driver increases the probability of survival
    if(new_geno['gnad'] >= 1) 
    {
      newval = psurv + delta_increase#/(M %>% filter(id == !!id) %>% pull(gnad))
      #if(newval >= 1) newval = 0.999
      return(newval)
    }
  }
  change_psurv = Vectorize(change_psurv, vectorize.args = 'psurv')
   
  psurv_parent = M %>% filter(id == !!parent_id) %>% pull(psurv)
   
  # True or false - is a new clone based on the psurv phenotype
  is_new_clone_surv = (new_geno['gnak'] >= 1 | new_geno['gnad'] >= 1) & (psurv_parent < 0.999 | psurv_parent > 0.001 )
  
  ######################################################################################################################################################
  ############################################# PHENOTYPE 2: Probability of immune attack ##############################################################
  
  
  #Function to modify the probability of immune attack
  change_pimmune = function(px)
  {
    # immunogenic mutations increase the chance to be detected by the IS
    if(new_geno['gnai'] >= 1)
    {
      newval2 = px + rgamma(1, 0.5, rate = 10)
      if(newval2 > 1) newval2 = 1
      return(newval2)
    }
    # escape mutations decrease the probability to be detected by the IS
    if(new_geno['gnae'] >= 1)
    {
      newval2 = 0 #- rgamma(1, 0.5, rate = 10)
      if(newval2 < 0) newval2 = 0
      return(newval2)
    }
    else{
    return(NULL)
    }
  }
  change_pimmune = Vectorize(change_pimmune, vectorize.args = 'px')
  
  pimmune_parent = M %>% filter(id == !!parent_id) %>% pull(pimmune)
  
  # True or false - is a new clone based on the pimmune phenotype
  is_new_clone_immune = (new_geno['gnai'] == 1) & pimmune_parent == 0
  
  #modify the phenotype
  M = M %>% 
    mutate(
      #survival probability
      psurv           = ifelse(id == !!id & is_new_clone_surv, change_psurv(psurv), psurv),
      #immune attack probability
      pimmune         = ifelse(id == !!id & is_new_clone_immune, change_pimmune(pimmune), pimmune)
    )
  
  
  #is_truly_newclone = is_new_clone_surv #| is_new_clone_immune
  
   ###set new strategy, new clone id and update parent clone id
  M = M %>%
    mutate(
      ##Strategy should always be updated
      #strategy        = ifelse(id == !!id & is_truly_newclone, outcome, strategy),
      strategy        = ifelse(id == !!id, outcome, strategy),
      clone_id        = ifelse(id == !!id & is_new_clone_surv, newclone_id, clone_id),
      clone_I_id        = ifelse(id == !!id & is_new_clone_immune, newclone_I_id, clone_I_id),
      parent_clone_id = ifelse(id == !!id & is_new_clone_surv, old_parent_clone_id, parent_clone_id)
    )
  M
}

###Get fitness increase 
#Gompertz function for increase
gompertz = function(a,b,c,t){
  val = a * exp(-b * (exp(-c*t)))
  return(val)
}

##Set a genotype for a cell
set_geno = function(M, id, g) {
  M = M %>%
    mutate(
      gnad = ifelse(id == !!id,  g['gnad'], gnad),
      gnap = ifelse(id == !!id,  g['gnap'], gnap),
      gnai = ifelse(id == !!id,  g['gnai'], gnai),
      gnae = ifelse(id == !!id,  g['gnae'], gnae),
      gnak = ifelse(id == !!id,  g['gnak'], gnak),
      gns  = ifelse(id == !!id,  g['gns'] , gns) ,
      gnsi = ifelse(id == !!id,  g['gnsi'], gnsi),
      gnsd = ifelse(id == !!id,  g['gnsd'], gnsd)
    )
  
  M
}

### 3 models for immunogenic attack
##Kill a cell given a probability of immune attack
set_immunedead1 = function(M, id, t, patak) {
  M %>%
    mutate(alive = ifelse(id == !!id & gnai > 0 & patak > runif(1) & gnae == 0,  FALSE, alive),
           t_e = ifelse(id == !!id,  t, t_e))
}

##Kill a cell with an increasing probability given a probability of immune attack
set_immunedead2 = function(M, id, t, patak) {
  M %>%
    mutate(alive = ifelse(id == !!id & pimmune > runif(1) & patak > runif(1) & gnae == 0,  FALSE, alive),
           t_e = ifelse(id == !!id,  t, t_e))
}

####Kill an immunogenic clone given a probability of immune attack
set_immunedead3 = function(M, t, pattack) {
  
  immunogenic_clones = immunogenic_cloneid(M, !!t)
  table_immunogenic = count_immunogenic_cells(M , !!t)
  
  if (length(immunogenic_clones) > 0 & (!!pattack > runif(1))) {
    icells=table_immunogenic %>% filter(clone_I_id %in% immunogenic_clones) %>% pull(count)
    print(paste(icells,
                paste(
                  paste("killed in Immunogenic clone",immunogenic_clones,sep=" "), paste("at time", !!t, sep=" ")
                )
                ,sep=" "))
    M = M %>%
      mutate(alive = ifelse(
        #(clone_id %in% immunogenic_clones) & (gnai > 0) & (patak > runif(1)),
        (clone_I_id %in% immunogenic_clones) & (gnai > 0) & (gnae == 0),
        FALSE,
        alive
    )
      )
  }else{
    M = M
   }
return(M)
}

##Flag a cell death and set the time of exit pool
set_dead_cell = function(M, id, t) {
  M %>%
    mutate(alive = ifelse(id == !!id,  FALSE, alive),
           t_e = ifelse(id == !!id,  t, t_e))
}

# Kill a cell
kill_cell = function(M, id, t) {
  set_dead_cell(
    M %>%
      mutate(
      strategy = ifelse(id == !!id,  "Death", strategy)
      ),
    id, 
    t
    )
}

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Computations functions that evaluate values
# from the state M of the model
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#select the configuration of mutations of a single cell when finish dividing
chooseStrategy <- function(M, id) {
  #Define outcome as Unknown first
  outcome = "Initial"
  #Accumulation of escape mutations
  if ((M %>% filter(id == !!id) %>% pull(gnae)) >= 1) {
    outcome = "Escape"
  }  
  else if ( (M %>% filter(id == !!id) %>% pull(gnai)) >= 1) {
    outcome = "Immunogenic"
  }
  ##Deleterious mutation, Deleterious mutation increases chance of dying
  #else if (M %>% filter(id == !!id) %>% pull(gnak) >= 1) {
  #  outcome = "Deleterious"
  #}
  else if ( (M %>% filter(id == !!id) %>% pull(gnad)) >= 1) {
    outcome = "Driver"
  }
  else{
    outcome = "Survive"
  }
  outcome
}

##Compute dnds for a sigle cell
compute_dnds = function(M, id) {
  g = get_geno(M, id)
  na = sum(g[1:5])
  ns = sum(g[6:8])
  dnds = na / (3 * ns)
  return(dnds)
}

##Compute dnds for a sigle cell (used for global values)
compute_gdnds = function(gnad, gnap, gnai, gnae, gnak, gns, gnsi, gnsd){
  na = gnad + gnap + gnai + gnae + gnak
  ns = gns + gnsi + gnsd
  dnds = na / (3 * ns)
  return(dnds)
}

##Compute dnds immunopeptidome for a single cell
compute_dnds_i = function(M, id) {
  g = get_geno(M, id)
  na = sum(g['gnai'])
  ns = sum(g['gnsi'])
  dnds_i = na / (3 * ns)
  return(dnds_i)
}

##Compute dnds immunopeptidom for global values)
compute_gdnds_i = function(gnai, gnsi) {
  na = gnai
  ns = gnsi
  dnds_i = na / (3 * ns)
  return(dnds_i)
}

##Compute dnds immunopeptidome for a single cell
compute_dnds_d = function(M, id) {
  g = get_geno(M, id)
  na = sum(g['gnad'])
  ns = sum(g['gnsd'])
  dnds_d = na / (3 * ns)
  return(dnds_d)
}

##Compute dnds immunopeptidom for global values)
compute_gdnds_d = function(gnad, gnsd) {
  na = gnad
  ns = gnsd
  dnds_d = na / (3 * ns)
  return(dnds_d)
}

##Sample a set of new somatic mutations for a cell
chooseMuts <- function(M, id, a0, N) {
  pc <- get_prob(M, id)
  ##Get a set of mutations based on the poisson of the expected
  Muts <- sample(a0, rpois(1, N), prob = pc, replace = T)
  
  g0 = get_geno(M, id)
  g0[T] = 0
  
  g1 <- table(Muts)
  g0[names(g1)] = g1
  g0
}

#Create a new cell from a parent cell
create_cell = function(M, new_id, parent_id) {
  f = M %>% filter(id == !!parent_id)
  f$id = new_id
  f$parent = parent_id
  f$alive = TRUE
  f$t_0 = f$t_0 + 1
  
  # Update clone info (you belong to the progeny of your ancestor's clone)
  bind_rows(M, f)
}

# Simulate a daugher cell at time t
simulate_daughter_cell = function(M, parent_id, t, a0, N) {
  # 1) new cell ID
  new_cell_id = get_new_id(M)
  
  # 2) Copy cell from the father
  M = create_cell(M, new_cell_id, parent_id)
  
  # 3) New genotype for this cell (acquired mutations + old mutations)
  # 3.1) choose mutations from probability vector
  new_geno = chooseMuts(M, new_cell_id, a0, N)
  # 3.2)  get current genotype and sum with old
  g_somatic = get_geno(M, parent_id) + new_geno
  # 3.3) Update genotype in memory
  M = set_geno(M, new_cell_id, g_somatic)
  
  # 4) Take new cell id, new genotype, and play the strategy
  M = playStrategy(M, new_cell_id, new_geno, parent_id)
  
  # 5) Update dnds, dnds_d and dndsi in memory and t_e
  M = M %>% 
    mutate(
      dnds = ifelse(id == !!new_cell_id,  compute_dnds(M, new_cell_id), dnds),
      dnds_d = ifelse(id == !!new_cell_id,  compute_dnds_d(M, new_cell_id), dnds_d),
      dnds_i = ifelse(id == !!new_cell_id,  compute_dnds_i(M, new_cell_id), dnds_i),
      t_e = ifelse(id == !!new_cell_id,  t, t_e)
  )
  
  # 6) Obtain new strategy and if new strategy 
  #new_strategy = M %>% filter(id == !!new_cell_id) %>% pull(strategy)

  #if (new_strategy == "Cell dies")
  #  M = set_dead_cell(M, new_cell_id, t + 1)
  #  M
}

count_immunogenic_cells = function(M, t) {
  M %>% 
    filter(
      t_0 == !!t,
      #pimmune >= 0.2
      strategy == "Immunogenic",
      gnae == 0,
      gnai > 0
      # alive
    ) %>%
    group_by(clone_I_id) %>%
    summarise(count=n())
}

immunogenic_cloneid = function(M, t) {
  M %>% 
    filter(
      t_0 == !!t,
      #pimmune >= 0.2
      strategy == "Immunogenic",
      gnae == 0,
      gnai > 0
      # alive
      ) %>%
    group_by(clone_I_id) %>%
    filter(n() > 50) %>%
    pull(clone_I_id) %>%
    unique()
}

####Functions for new dnds calculation
get_t<-function(M,id){
  time<-M %>% filter(id==!!id) %>% pull(t_0)
  return(time)
}

###Obtain the daugther ID
get_daughter_id<-function(M,id){
  t_parent<-get_t(M,id)
  t_child<-t_parent+1
  childid<-M %>% filter(t_0==t_child,parent==!!id,strategy!="Death") %>% pull(id)
  return(childid)
}

##Obtain a random ID
getrandomid<-function(x){
  mutid<-runif(x,min=0 ,max=50000000)
  return(mutid)
}

##Assign new mutation ID
assign_newmutids<-function(vector){
  vector2<-tibble::enframe(vector)
  vector2<-vector2 %>% mutate(newid=ifelse(value==TRUE,paste(name,as.integer(getrandomid(1)),sep="_"),"NA")) %>% select(newid) %>% filter(newid!="NA") %>% pull(newid)
  return(vector2)
}

##Compare genotype of child versus paternal id
compare_genotypes<-function(paternal,child){
  oldids=""
  for (i in 0:10){
    vector<- (child-paternal > i)
    newids<-assign_newmutids(vector)
    oldids=c(oldids,newids)
  }
  oldids<-oldids[oldids != ""]
  oldids<-unique(oldids)
  return(oldids)
}

## Obtain full genotype for each cell in the model
get_fullgenotype<-function(i,MHF){
  fullgenotype <- paste(MHF %>% filter(cellid %in% i) %>% pull(Mut_IDs),collapse = ":") %>% str_replace_all(.,":+",":")
  return(fullgenotype)
}

#####Main function to obtain sequencing table
get_sequencing<-function(x){
  
  M<-x
  
  MHF<-tibble( t_0 =NA, cellid = NA, paternalID = NA, Mut_IDs = NA)
  
  for(id in 1:length(M$id)){
    CHILD<-tibble(t_0 = NA, cellid = NA, paternalID = NA, Mut_IDs = NA)
    t_actual<-get_t(M,id)
    c_actual<-id
    
    #Get paternal genotype
    paternal_geno<-get_geno(M,id)
    
    #Get child genotypes
    childs<-get_daughter_id(M,id)
    
    #Populate new matrix with the information of the daugther and the new mut ids
    if(length(childs)>0){
      for (j in childs){
        child_geno<-get_geno(M,j)
        newmutsids<-compare_genotypes(paternal = paternal_geno,child = child_geno)
        CHILD<-tibble(t_0 = t_actual+1, cellid = j, paternalID =  c_actual, Mut_IDs = paste(newmutsids,collapse = ":"))
        MHF = bind_rows(MHF,CHILD)
      }
    }
  }
  

  
  ####### Pull table with mutation ids, number of cells, frequency, and type
  # Call id the "to" and parent the "from" of each edge
  M = M %>%
    mutate(to = id, from = parent)
  
  #Example
  cells = M %>% filter(t_0==max(t_0),strategy!="Death") %>% pull(to)
  tfinal <- M %>%  filter(t_0==max(t_0)) %>% pull(t_0) %>% unique(.)
  total  <- length(cells)
  
  if(total == 0){
    print("All cells dead at last time point")
    final_seq_unfiltered<-NA
    return(final_seq_unfiltered)
    next
  }else{
    print("Cells are alive in last generation, compute dN/dS")
  }
  #Get values for each cell
  Mad = ctree:::DataFrameToMatrix(M)
  # 2) Transpose it so that now every edge goes from
  # child to father
  t_Mad = t(Mad)
  # 3)
  t_M = ctree:::MatrixToDataFrame(t_Mad)
  ancestors<-lapply(cells, ctree:::reach, df = t_M)
  
  sequencing<-""
  for (a in ancestors){
    temp_genotype<-get_fullgenotype(a,MHF)
    sequencing<-paste(sequencing,temp_genotype,sep=":") %>% str_replace_all(.,":+",":")
  }
  
  sequencing_vec<-unlist(strsplit(sequencing, split=":"))
  
  final_seq_unfiltered<-as_tibble(table(sequencing_vec)) %>% mutate(frequency=n/total)
  final_seq_unfiltered<-final_seq_unfiltered %>% mutate(muttype=(ifelse(str_detect(sequencing_vec, 'gna'),"Non-syn","syn")))
  
  return(final_seq_unfiltered)
}

## Get frequency dnds
calc_freq_dnds<-function(x,sequencing,freq=0.01){
  summary_dnds_table=tibble(simulation=NA,pop=NA,time=NA,mut_na=NA,mut_ns=NA,mut_nad=NA,
                            mut_nsd=NA,mut_nai=NA,mut_nsi=NA,mut_nae=NA,dnds_global=NA,
                            dnds_driver=NA,dnds_immune=NA,min_freq=NA)
  
  for (i in 1:length(x)){
    #print(paste("Simulation_nr",i,sep ="_"))
    #Iterate through simulations
    M <- x[[i]]

    if(unique(is.na(sequencing[[i]][1]))){ 
      #print(paste("No Sequencing data in simulation ",i,sep =""))
      next}
    #Final table data
    temp_dnds_table=tibble(simulation=NA,pop=NA,time=NA,mut_na=NA,mut_ns=NA,mut_nad=NA,
                           mut_nsd=NA,mut_nai=NA,mut_nsi=NA,mut_nae=NA,dnds_global=NA,
                           dnds_driver=NA,dnds_immune=NA,min_freq=!!freq)
    
    #Assign simulation number
    temp_dnds_table <- temp_dnds_table %>% mutate(simulation = i)
    
    #Define cells to query, time to query, and total population size to query
    cells = M %>% filter(t_0==max(t_0),strategy!="Death") %>% pull(id)
    tfinal <- M %>%  filter(t_0==max(t_0)) %>% pull(t_0) %>% unique(.)
    total  <- length(cells)
    
    #If population is extinct do not compute dN/dS
    if(total == 0){next}
    
    #Get matrix of cells and explicit genotypes
    final_seq_unfiltered<-dplyr::as_tibble(sequencing[[i]])
    
    final_seq<- final_seq_unfiltered %>% filter(frequency>=!!freq)
    
    #get_dnds_global
    na<-as.integer(count(final_seq %>% rename(n_cells=n) %>% filter(str_detect(sequencing_vec, 'gna'))))
    ns<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gns'))))
    dnds<-na/(3*ns)
    
    #get_dnds_driver
    nad<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnad'))))
    nsd<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnsd'))))
    dndsD<-nad/(3*nsd)
    
    #get_dnds_Immune
    nai<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnai'))))
    nsi<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnsi'))))
    dndsI<-nai/(3*nsi)
    
    nae<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnae'))))
    
    temp_dnds_table = temp_dnds_table %>% mutate(pop = total,
                                                 time=tfinal,
                                                 mut_na=na,mut_ns=ns,
                                                 mut_nad=nad,mut_nsd=nsd,
                                                 mut_nai=nai,mut_nsi=nsi,
                                                 mut_nae=nae,
                                                 dnds_global=dnds,
                                                 dnds_driver=dndsD,
                                                 dnds_immune=dndsI,
                                                 min_freq=freq )
    
    
    summary_dnds_table = bind_rows(summary_dnds_table,temp_dnds_table)
    
  }
  summary_dnds_table<-summary_dnds_table[complete.cases(summary_dnds_table$simulation),]
  return(summary_dnds_table)
}

## Calculate if escape is clonal or subclonal
calc_clonal_escape<-function(x,sequencing,freq=0.01,escape_freq_min=0.1,escape_freq_max=1){
  summary_dnds_table=tibble(simulation=NA,pop=NA,time=NA,mut_na=NA,mut_ns=NA,mut_nad=NA,
                            mut_nsd=NA,mut_nai=NA,mut_nsi=NA,mut_nae=NA,dnds_global=NA,
                            dnds_driver=NA,dnds_immune=NA,min_freq=NA)
  
  for (i in 1:length(x)){
    #print(paste("Simulation_nr",i,sep ="_"))
    #Iterate through simulations
    M <- x[[i]]
    
    if(unique(is.na(sequencing[[i]][1]))){ 
      #print(paste("No Sequencing data in simulation ",i,sep =""))
      next}
    #Final table data
    temp_dnds_table=tibble(simulation=NA,pop=NA,time=NA,mut_na=NA,mut_ns=NA,mut_nad=NA,
                           mut_nsd=NA,mut_nai=NA,mut_nsi=NA,mut_nae=NA,dnds_global=NA,
                           dnds_driver=NA,dnds_immune=NA,min_freq=!!freq)
    
    #Assign simulation number
    temp_dnds_table <- temp_dnds_table %>% mutate(simulation = i)
    
    #Define cells to query, time to query, and total population size to query
    cells = M %>% filter(t_0==max(t_0),strategy!="Death") %>% pull(id)
    tfinal <- M %>%  filter(t_0==max(t_0)) %>% pull(t_0) %>% unique(.)
    total  <- length(cells)
    
    #If population is extinct do not compute dN/dS
    if(total == 0){next}
    
    #Get matrix of cells and explicit genotypes
    final_seq_unfiltered<-sequencing[[i]]
    
    final_seq<- dplyr::as_tibble(final_seq_unfiltered %>% filter(frequency>=!!freq))
  
    max_freq_escape<- unique(final_seq_unfiltered %>% filter(.,str_detect(sequencing_vec,"gnae")) %>% filter(frequency==max(frequency)) %>% pull(frequency))
  
    if(purrr::is_empty(max_freq_escape)){next}
    
    if(max_freq_escape>=escape_freq_min & max_freq_escape<escape_freq_max){
      #get_dnds_global
      na<-as.integer(count(final_seq %>% rename(n_cells=n) %>% filter(str_detect(sequencing_vec, 'gna'))))
      ns<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gns'))))
      dnds<-na/(3*ns)
      
      #get_dnds_driver
      nad<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnad'))))
      nsd<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnsd'))))
      dndsD<-nad/(3*nsd)
      
      #get_dnds_Immune
      nai<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnai'))))
      nsi<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnsi'))))
      dndsI<-nai/(3*nsi)
      
      nae<-as.integer(count(final_seq %>% rename(n_cells=n) %>%filter(str_detect(sequencing_vec, 'gnae'))))
      temp_dnds_table = temp_dnds_table %>% mutate(pop = total,
                                                   time=tfinal,
                                                   mut_na=na,mut_ns=ns,
                                                   mut_nad=nad,mut_nsd=nsd,
                                                   mut_nai=nai,mut_nsi=nsi,
                                                   mut_nae=nae,
                                                   dnds_global=dnds,
                                                   dnds_driver=dndsD,
                                                   dnds_immune=dndsI,
                                                   min_freq=freq )
      
      
      summary_dnds_table = bind_rows(summary_dnds_table,temp_dnds_table)
    }else{next}
  }
  summary_dnds_table<-summary_dnds_table[complete.cases(summary_dnds_table[,1]),]
  return(summary_dnds_table)
}

## Estimate confidence interval for cohort analysis
calc_freq_ci<-function(x){
  
  sum_x <-
    x %>% dplyr::group_by(immune, escape, min_freq) %>% dplyr::summarise(
      sum_na  = sum(mut_na),
      sum_ns  = sum(mut_ns),
      sum_nad = sum(mut_nad),
      sum_nsd = sum(mut_nsd),
      sum_nai = sum(mut_nai),
      sum_nsi = sum(mut_nsi),
      sum_nae = sum(mut_nae),
      count = n(),
      dnds_global = (sum(mut_na) /
                       (3 * sum(mut_ns))),
      dnds_g_low = get_lowCI_freq(sum_na,3*sum_ns,count,count),
      dnds_g_high = get_highCI_freq(sum_na,3*sum_ns,count,count),
      
      dnds_driver = (sum(mut_nad) /
                       (3 * sum(mut_nsd))),
      dnds_d_low = get_lowCI_freq(sum_nad,3*sum_nsd,count,count),
      dnds_d_high = get_highCI_freq(sum_nad,3*sum_nsd,count,count),
      dnds_immune = (sum(mut_nai) /
                       (3 * sum(mut_nsi))),
      dnds_i_low = get_lowCI_freq(sum_nai,3*sum_nsi,count,count),
      dnds_i_high = get_highCI_freq(sum_nai,3*sum_nsi,count,count),
      freq = mean(min_freq)
    )
  
  return(sum_x)
  }

## Get low CI
get_lowCI_freq=function(x1,x2,n,m){
  library(dpcR)
  z<-rateratio.test::rateratio.test(c(x1,x2),c(n,m))
  return(z$conf.int[1])
}

##get High CI
get_highCI_freq=function(x1,x2,n,m){
  #library(dpcR)
  z<-rateratio.test::rateratio.test(c(x1,x2),c(n,m))
  return(z$conf.int[2])
}

##Set parameters for plotting cohort analysis
set_params<-function(x,mod,imm,esc){
  x <- x %>% mutate(model=mod,immune=imm,escape=esc)
}  

##New function to get clonal escape id and get runs with that data
get_clonal_escape_simID=function(dnds_table,x){
  summary_dnds_table=tibble(simulation=NA,pop=NA,time=NA,mut_na=NA,mut_ns=NA,mut_nad=NA,
                            mut_nsd=NA,mut_nai=NA,mut_nsi=NA,mut_nae=NA,dnds_global=NA,
                            dnds_driver=NA,dnds_immune=NA,min_freq=NA,model=NA,immune=NA,escape=NA)
  
  temp_dnds_table=NULL
  tmp_positive_clone = list()
  #output_list = list()
  for (i in 1:length(x)){
    
    if(unique(is.na(x[[i]][1]))){ 
      #print(paste("No Sequencing data in simulation ",i,sep =""))
      next}  
    
    
    tmp_positive_clone[[i]]=x[[i]] %>% dplyr::filter(str_detect(sequencing_vec,"gnae")) %>% dplyr::filter(frequency>0.5)
    
    if(dplyr::count(tmp_positive_clone[[i]])>0){
      ##Create temporary tibble
      temp_dnds_table = tibble(simulation=NA,pop=NA,time=NA,mut_na=NA,mut_ns=NA,mut_nad=NA,
                               mut_nsd=NA,mut_nai=NA,mut_nsi=NA,mut_nae=NA,dnds_global=NA,
                               dnds_driver=NA,dnds_immune=NA,min_freq=NA,model=NA,immune=NA,escape=NA)
      ##Fill temporary tibble
      temp_dnds_table = dnds_table %>% dplyr::filter(simulation == i)
      ##remove NA rows
      summary_dnds_table = bind_rows(summary_dnds_table,temp_dnds_table)
    }
    
    
  }
  summary_dnds_table <- summary_dnds_table %>% dplyr::filter(!is.na(simulation))
  return(summary_dnds_table)
}


