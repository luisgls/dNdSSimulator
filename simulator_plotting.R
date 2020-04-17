

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Plotting functions
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


###Different functions for summarising and plotting the results of the simulations
##plot multiple plots from single simulation and save
plotter = function(j, M)
{
  library(tidyverse)
  source("functions_ZCsimulator.R")
  
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
        add_clone_tree = TRUE
      )
    
    ggsave(
      p1,
      filename = paste(j, "simulation.pdf", sep = "_"),
      width = 35,
      height = 15,
      units = "cm"
    )
    
    p2 <-
      ggstatsplot::ggbetweenstats(
        data = M %>% filter(t_0 == max(t_0), dnds != Inf),
        x = strategy,
        y = dnds,
        pairwise.comparisons = TRUE,
        type = "np",
        p.adjust.method = "BH",
        # palette = "default_npg", # choosing a different color palette
        xlab = "",
        ylab = "dN/dS",
        # label for the x-axis variable
        messages = FALSE
      ) + geom_abline(slope = 0,
                      intercept = 1,
                      colour = "red")
    
    ggsave(
      p2,
      filename = paste(j, "dnds_simulation.pdf", sep = "_"),
      width = 15,
      height = 15,
      units = "cm"
    )
  }
}

##plot multiple plots from single simulation but dont save
plotter_nosave = function(M)
{
  library(tidyverse)
  source("functions_ZCsimulator.R")
  
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
        add_clone_tree = TRUE
      )
    p1
    # ggsave(
    #   p1,
    #   filename = paste(j, "simulation_d01_i01.pdf", sep = "_"),
    #   width = 35,
    #   height = 15,
    #   units = "cm"
    # )
    
    p2 <-
      ggstatsplot::ggbetweenstats(
        data = M %>% filter(t_0 == max(t_0), dnds != Inf),
        x = strategy,
        y = dnds,
        pairwise.comparisons = TRUE,
        type = "np",
        p.adjust.method = "BH",
        # palette = "default_npg", # choosing a different color palette
        xlab = "",
        ylab = "dN/dS",
        # label for the x-axis variable
        messages = FALSE
      ) + geom_abline(slope = 0,
                      intercept = 1,
                      colour = "red")
    
    # ggsave(
    #   p2,
    #   filename = paste(j, "dnds_simulation_d01_i01.pdf", sep = "_"),
    #   width = 15,
    #   height = 15,
    #   units = "cm"
    #)
    p2
  }
}

##Create a table statistic of dN/dS for different simulations
stats_dnds <- function(x)
{
  fulltimeseries <-
    tibble(
      t_0 = NA,
      median_dnds = NA,
      median_dndsi = NA,
      median_dndsd = NA,
      dnds = NA,
      dnds_lci = NA,
      dnds_hci = NA,
      dnds_i = NA,
      dnds_ilci = NA,
      dnds_ihci = NA,
      dnds_d = NA,
      dnds_dlci = NA,
      dnds_dhci = NA,
      count = NA,
      expon = NA,
      sum_na = NA,
      sum_ns = NA,
      gnad = NA,
      gnap = NA,
      gnai = NA,
      gnae = NA,
      gnak = NA,
      gns = NA,
      gnsi = NA,
      gnsd = NA,
      numclones = NA
    )
  fulltimeseries <- fulltimeseries[complete.cases(fulltimeseries), ]
  
  for (i in 1:length(x)) {
    M <- x[[i]]
    M = M %>% mutate(t_0 = ifelse(is.na(t_0), max(M$t_0) + 1, t_0))
    M = M %>% filter(strategy != "Death")
    
    pop = M %>%
      group_by(t_0) %>%
      summarise(count = n()) %>%
      mutate(expon = 2 ^ (t_0))
    
    ####Get the population parameters for the population dN/dS
    DNDS = M %>%
      group_by(t_0) %>%
      summarise(
        gnad = sum(gnad),
        gnap = sum(gnap),
        gnai = sum(gnai),
        gnae = sum(gnae),
        gnak = sum(gnak),
        gns = sum(gns),
        gnsi = sum(gnsi),
        gnsd = sum(gnsd),
        median_dnds = median(dnds),
        median_dndsi = median(dnds_i),
        median_dndsd = median(dnds_d)
      ) %>%
      mutate(
        dnds = compute_gdnds(gnad, gnap, gnai, gnae, gnak, gns, gnsi, gnsd),
        dnds_i = compute_gdnds_i(gnai, gnsi),
        dnds_d = compute_gdnds_d(gnad, gnsd),
        dnds_lci =  dnds - (1.96 * (sqrt(
          dnds / (gnad + gnap + gnai + gnae + gnak + gns + gnsi)
        ))),
        dnds_hci =  dnds + (1.96 * (sqrt(
          dnds / (gnad + gnap + gnai + gnae + gnak + gns + gnsi)
        ))),
        dnds_ilci = dnds_i - (1.96 * (sqrt(
          dnds_i / (gnai + gnsi)
        ))),
        dnds_ihci = dnds_i + (1.96 * (sqrt(
          dnds_i / (gnai + gnsi)
        ))),
        dnds_dlci = dnds_d - (1.96 * (sqrt(
          dnds_d / (gnad + gnsd)
        ))),
        dnds_dhci = dnds_d + (1.96 * (sqrt(
          dnds_d / (gnad + gnsd)
        )))
      )
    
    ##Get the totals
    totals = M %>%
      group_by(t_0) %>%
      summarise(
        gnad = sum(gnad),
        gnap = sum(gnap),
        gnai = sum(gnai),
        gnae = sum(gnae),
        gnak = sum(gnak),
        gns = sum(gns),
        gnsi = sum(gnsi),
        gnsd = sum(gnsd),
        sum_na = gnad + gnap + gnai + gnae + gnak,
        sum_ns = gns + gnsi + gnad,
        numclones = length(unique(clone_id))
      )
    
    ###In order to plot we need to join a simplified table
    time_series = full_join(pop, DNDS) %>% full_join(totals) %>%
      dplyr::select(
        t_0,
        median_dnds,
        median_dndsi,
        median_dndsd,
        dnds,
        dnds_lci,
        dnds_hci,
        dnds_i,
        dnds_ilci,
        dnds_ihci,
        dnds_d,
        dnds_dlci,
        dnds_dhci,
        count,
        expon,
        sum_na,
        sum_ns,
        gnad,
        gnap,
        gnai,
        gnae,
        gnak,
        gns,
        gnsi,
        gnsd,
        numclones
      )
    
    time_series = time_series %>% mutate(simulation = i)
    
    fulltimeseries <- bind_rows(fulltimeseries, time_series)
    
  }
  
  return(fulltimeseries)
}

##Obtain summary statistics table for Expansion, extinction and ongoing
evodevo = function(x, min, max)
{
  fullevolution <-
    tibble(
      simulation_number = NA,
      time = NA,
      event = NA,
      fin_dnds = NA,
      fin_dndsd = NA,
      fin_dndsi = NA
    )
  #fullevolution <- fullevolution[complete.cases(fullevolution), ]
  
  for (i in 1:length(x)) {
    M <- x[[i]]
    M = M %>% mutate(t_0 = ifelse(is.na(t_0), max(M$t_0) + 1, t_0))
    #M = M %>% filter(strategy != "Death")
    popsize <- nrow(M %>% filter(alive, t_0 == max(M %>% dplyr::select(t_0))))
    lasttime <- max(M %>% dplyr::select(t_0))
    finaldnds <-
      M %>% filter(t_0 == (lasttime)) %>% summarise(finaldnds = (sum(gnad) +
                                                                   sum(gnap) + sum(gnai) + sum(gnae) + sum(gnak)) / (3 * (sum(gns) + sum(gnsi) +
                                                                                                                            sum(gnsd)))) %>% pull(finaldnds)
    finaldndsd <-
      M %>% filter(t_0 == (lasttime)) %>% summarise(finaldndsd = (sum(gnad)) / (3 *
                                                                                  (sum(gnsd)))) %>% pull(finaldndsd)
    finaldndsi <-
      M %>% filter(t_0 == (lasttime)) %>% summarise(finaldndsi = (sum(gnai)) / (3 *
                                                                                  (sum(gnsi)))) %>% pull(finaldndsi)
    counts <-
      M %>% filter(t_0 == (lasttime)) %>% summarise(counts = sum(gnad)) %>% pull (counts)
    counts2 <-
      M %>% filter(t_0 == (lasttime)) %>% summarise(counts2 = sum(gnai)) %>% pull (counts2)
    
    finaldndsd = ifelse(is.na(finaldndsd), 0, finaldndsd)
    finaldndsd = ifelse(is.infinite(finaldndsd), counts, finaldndsd)
    finaldndsi = ifelse(is.na(finaldndsi), 0, finaldndsi)
    finaldndsi = ifelse(is.infinite(finaldndsi), counts2, finaldndsi)
    
    evolution <- tibble(
      simulation_number = i,
      time = lasttime,
      event = ifelse(
        popsize > max,
        "Expansion",
        ifelse(popsize < min, "Extinction", "Ongoing")
      ),
      fin_dnds = finaldnds,
      fin_dndsd = finaldndsd,
      fin_dndsi = finaldndsi
    )
    fullevolution = bind_rows(fullevolution, evolution)
    
  }
  fullevolution <- fullevolution[complete.cases(fullevolution), ]
  return(fullevolution)
}

##summarise final dnds for different simulations
dnds_time = function(x, min, max)
{
  fullevolution <- tibble(simulation_number = NA,
                          time = NA,
                          event = NA)
  fullevolution <- fullevolution[complete.cases(fullevolution), ]
  
  for (i in 1:length(x)) {
    M <- x[[i]]
    M = M %>% dplyr::mutate(t_0 = ifelse(is.na(t_0), max(M$t_0) + 1, t_0))
    M = M %>% filter(strategy != "Death")
    popsize <- nrow(M %>% filter(alive, t_0 == max(M %>% dplyr::select(t_0))))
    lasttime <- max(M %>% dplyr::select(t_0))
    finaldnds <-
      M %>% group_by(t_0) %>% summarise(
        simulation_number = i,
        finaldnds = (sum(gnad) + sum(gnap) + sum(gnai) + sum(gnae) + sum(gnak)) / (3 *
                                                                                     (sum(gns) + sum(gnsi) + sum(gnsd)))
      ) %>% pull(finaldnds)
    finaldndsd <-
      M %>% group_by(t_0) %>% summarise(simulation_number = i,
                                        finaldndsd = (sum(gnad)) / (3 * (sum(gnsd))))
    
    
    evolution <- tibble(
      simulation_number = i,
      time = lasttime,
      event = ifelse(
        popsize > max,
        "Expansion",
        ifelse(popsize < min, "Extinction", "Ongoing")
      )
    )
    
    evojoin <- full_join(evolution, finaldnds, by = simulation_number)
    fullevolution = bind_rows(fullevolution, evolution)
    
    
  }
  return(fullevolution)
}

##plot multiple plots from single simulation and save
plotter2 = function(x, name)
{
  fulltimeseries = x
  
  p1 <- ggstatsplot::ggbetweenstats(
    data = fulltimeseries %>% filter(sum_ns > 0),
    x = t_e,
    y = dnds,
    #grouping.var = immune,
    #pairwise.comparisons = TRUE,      # display significant pairwise comparisons
    #pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
    type = "np",
    p.adjust.method = "BH",
    # method for adjusting p-values for multiple comparisons
    #bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
    #package = "ggsci", # package from which color palette is to be taken
    #palette = "nrc_npg", # choosing a different color palette
    xlab = "",
    # label for the x-axis variable
    messages = FALSE
  ) + geom_abline(slope = 0,
                  intercept = 1,
                  colour = "red") #+  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(colourCount))
  
  p12 <- ggstatsplot::ggbetweenstats(
    data = fulltimeseries %>% filter(gnsi > 0),
    x = t_e,
    y = dnds_i,
    #grouping.var = immune,
    #pairwise.comparisons = TRUE,      # display significant pairwise comparisons
    #pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
    type = "np",
    p.adjust.method = "BH",
    # method for adjusting p-values for multiple comparisons
    #bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
    package = "ggsci",
    # package from which color palette is to be taken
    palette = "nrc_npg",
    # choosing a different color palette
    xlab = "",
    # label for the x-axis variable
    messages = FALSE
  ) + geom_abline(slope = 0,
                  intercept = 1,
                  colour = "red")
  
  p13 <- ggstatsplot::ggbetweenstats(
    data = fulltimeseries %>% filter(gnsd > 0),
    x = t_e,
    y = dnds_d,
    #grouping.var = immune,
    #pairwise.comparisons = TRUE,      # display significant pairwise comparisons
    #pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
    type = "np",
    p.adjust.method = "BH",
    # method for adjusting p-values for multiple comparisons
    #bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
    package = "ggsci",
    # package from which color palette is to be taken
    palette = "nrc_npg",
    # choosing a different color palette
    xlab = "",
    # label for the x-axis variable
    messages = FALSE
  ) + geom_abline(slope = 0,
                  intercept = 1,
                  colour = "red")
  
  p2 <- ggstatsplot::ggbetweenstats(
    data = fulltimeseries %>% filter(sum_ns > 0),
    x = t_e,
    y = sum_ns,
    #grouping.var = immune,
    #pairwise.comparisons = TRUE,      # display significant pairwise comparisons
    #pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
    type = "np",
    p.adjust.method = "BH",
    # method for adjusting p-values for multiple comparisons
    #bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
    package = "ggsci",
    # package from which color palette is to be taken
    palette = "nrc_npg",
    # choosing a different color palette
    #package = "ggsci",
    xlab = "",
    # label for the x-axis variable
    messages = FALSE
  ) #+  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(colourCount))
  
  p3 <- ggstatsplot::ggbetweenstats(
    data = fulltimeseries %>% filter(sum_ns > 0),
    x = t_e,
    y = count,
    #grouping.var = immune,
    #pairwise.comparisons = TRUE,      # display significant pairwise comparisons
    #pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
    type = "np",
    p.adjust.method = "BH",
    # method for adjusting p-values for multiple comparisons
    #bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
    #conf.level = 0.95,                # changing confidence level to 99%,
    #ggtheme = ggthemes::theme_fivethirtyeight(), # choosing a different theme
    #ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
    package = "ggsci",
    # package from which color palette is to be taken
    palette = "nrc_npg",
    # choosing a different color palette
    xlab = "Generation",
    # label for the x-axis variable
    messages = FALSE
  ) #+  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(colourCount))
  
  p4 <- ggstatsplot::ggbetweenstats(
    data = fulltimeseries %>% filter(sum_ns > 0) %>% mutate(norm = log(dnds) /
                                                              (count)),
    x = t_e,
    y = norm,
    #grouping.var = immune,
    #pairwise.comparisons = TRUE,      # display significant pairwise comparisons
    #pairwise.annotation = "p.value",  # how do you want to annotate the pairwise comparisons
    type = "np",
    p.adjust.method = "BH",
    # method for adjusting p-values for multiple comparisons
    #bf.message = TRUE,                # display Bayes Factor in favor of the null hypothesis
    #conf.level = 0.95,                # changing confidence level to 99%,
    #ggtheme = ggthemes::theme_fivethirtyeight(), # choosing a different theme
    #ggstatsplot.layer = FALSE, # turn off ggstatsplot theme layer
    package = "ggsci",
    # package from which color palette is to be taken
    palette = "nrc_npg",
    # choosing a different color palette
    xlab = "Generation",
    # label for the x-axis variable
    messages = FALSE
  )
  
  finalp1 <- ggarrange(p1, p12, p13, p2, p3)
  finalp1
  ggsave(
    paste("simul_", name, ".pdf" , sep = ""),
    scale = 4,
    width = 30,
    height = 25,
    units = "cm"
  )
  
  plot(log(fulltimeseries$dnds_d), log10(fulltimeseries$count))
}

##Plot global dnds, driver dnds and imnune dnds over time for single simulation
plot_five = function (x)
{
  medians <- aggregate(dnds ~  t_0, x, median)
  medians <- round(medians, 2)
  p1 <- ggboxplot(
    x,
    x = "t_0",
    y = "dnds",
    color = "#4DBBD5FF",
    palette = "jco"
  ) +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") +
    labs(title = "Overall selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal() + geom_text_repel(
      data = subset(medians, dnds < 1)
      ,
      aes(label = dnds, x = t_0 - 1),
      nudge_y = -0.5,
      nudge_x = 0.5,
      direction = "y"
    ) +
    geom_text_repel(
      data = subset(medians, dnds >= 1)
      ,
      aes(label = dnds, x = t_0 - 1),
      nudge_y = 0.5,
      nudge_x = 0.5,
      direction = "x"
    )
  
  medians <- aggregate(dnds_d ~  t_0, x, median)
  medians <- round(medians, 2)
  p2 <- ggboxplot(
    x,
    x = "t_0",
    y = "dnds_d",
    color = "#4DBBD5FF",
    palette = "jco"
  ) +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") + labs(title = "Driver selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal() + geom_text_repel(
      data = subset(medians, dnds_d < 1)
      ,
      aes(label = dnds_d, x = t_0 - 1),
      nudge_y = -0.5,
      nudge_x = 0.5,
      direction = "y"
    ) +
    geom_text_repel(
      data = subset(medians, dnds_d >= 1)
      ,
      aes(label = dnds_d, x = t_0 - 1),
      nudge_y = 0.5,
      nudge_x = 0.5,
      direction = "x"
    )
  
  medians <- aggregate(dnds_i ~  t_0, x, median)
  medians <- round(medians, 2)
  p3 <- ggboxplot(
    x,
    x = "t_0",
    y = "dnds_i",
    color = "#4DBBD5FF",
    palette = "jco"
  ) +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") + labs(title = "Immune selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal() + geom_text_repel(
      data = subset(medians, dnds_i < 1)
      ,
      aes(label = dnds_i, x = t_0 - 1),
      nudge_y = -0.5,
      nudge_x = 0.5,
      direction = "y"
    ) +
    geom_text_repel(
      data = subset(medians, dnds_i >= 1)
      ,
      aes(label = dnds_i, x = t_0 - 1),
      nudge_y = 0.5,
      nudge_x = 0.5,
      direction = "x"
    )
  
  p4 <- ggboxplot(
    x,
    x = "t_0",
    y = "count",
    color = "#00A087FF",
    palette = "jco"
  ) +
    geom_abline(slope = 0,
                intercept = 32,
                col = "red") +  labs(title = "Population size", y = 'Pop size', x = "Generation") +
    theme_minimal()
  
  p5 <- ggboxplot(x,
                  x = "t_0",
                  y = "gnap + gns",
                  color = "#3C5488FF")  +  labs(title = "Mutations", y = 'Muts', x = "Generation") +
    theme_minimal()
  
  finalp1 <- ggarrange(p4, p1, p2, p3, p5, nrow = 5, labels = "AUTO")
  finalp1
}

##Plot boxplots for dnds, pop size and mutation number in a boxplot
plot_three = function (x)
{
  medians <- aggregate(dnds ~  t_0, x, median)
  medians <- round(medians, 2)
  p1 <- ggboxplot(
    x,
    x = "t_0",
    y = "dnds",
    color = "#4DBBD5FF",
    palette = "jco"
  ) +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") +
    labs(title = "Overall selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal() + #+ stat_summary(fun.y=round(median,2), colour="darkred", geom="point",
    #   shape=1, size=3) +
    geom_text_repel(
      data = subset(medians, dnds < 1)
      ,
      aes(label = dnds, x = t_0),
      nudge_y = Inf,
      #nudge_x = 0.2,
      direction = "y",
      segment.size = 0.1,
      segment.alpha = 0.5
    ) +
    geom_text_repel(
      data = subset(medians, dnds >= 1)
      ,
      aes(label = dnds, x = t_0),
      nudge_y = Inf,
      #nudge_x = 0.2,
      direction = "x",
      segment.size = 0.1,
      segment.alpha = 0.5
    )
  p4 <- ggboxplot(
    x,
    x = "t_0",
    y = "count",
    color = "#00A087FF",
    palette = "jco"
  ) +
    geom_abline(slope = 0,
                intercept = 32,
                col = "red") +
    labs(title = "Population size", y = 'Pop size', x = "Generation") +
    theme_minimal()
  
  p5 <- ggboxplot(x,
                  x = "t_0",
                  y = "gnap + gns",
                  color = "#3C5488FF")  +
    labs(title = "Mutations", y = 'Muts', x = "Generation") +
    theme_minimal()
  
  finalp1 <- ggarrange(p4, p1, p5, nrow = 3, labels = "AUTO")
  finalp1
}

##Plot pop size, mutation number and overall dnds over time for multiple simulations
plot_threeV2 = function (x,
                         subset_n = 20,
                         maxt = 100,
                         maintitle = "Summary simulation")
{
  listofsimulations <- sample(1:max(x$simulation), subset_n)
  x <- x %>% filter(simulation %in% listofsimulations)
  
  
  medians <- aggregate(dnds ~  t_0, x, median)
  medians <- round(medians, 2)
  
  x$simulation <- as.factor(x$simulation)
  
  p1 <-
    ggplot(x, aes(
      x = t_0,
      y = dnds,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    geom_abline(slope = 0,
                intercept = 0,
                col = "red") +
    labs(title = "Overall selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") + scale_y_log10() + xlim(0, maxt)
  
  p2 <-
    ggplot(x, aes(
      x = t_0,
      y = count,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    labs(title = "Population size", y = 'Pop size', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") + scale_y_log10() +  geom_abline(slope = 0,
                                                                                       intercept = log10(16),
                                                                                       col = "red") + xlim(0, maxt)
  
  p3 <-
    ggplot(x, aes(
      x = t_0,
      y = sum_na,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    labs(title = "Mutations", y = 'Non-silent muts', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") + scale_y_log10() + xlim(0, maxt)
  
  finalp1 <- ggarrange(p2, p3, p1, nrow = 3, labels = "AUTO")
  
  annotate_figure(finalp1, top = text_grob(
    paste("Fig", maintitle, sep = " "),
    color = "black",
    face = "bold",
    size = 14
  ))
  
}

##Plot ongoing versus collapsse versus extinction for multiple simulations
plot_threeV3 = function (fulltimeevo,
                         subset_n = 50,
                         maxt = 100,
                         maintitle = "Summary simulation")
{
  #p002<-ggplot(data=fulltimeevo, aes(fill=event, x = time)) + geom_histogram(aes(alpha=0.5)) + scale_alpha(guide = 'none')  + theme_classic2()
  #ggviolin(fulltimeevo, x = "event", y = "time", fill = "event",
  #         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  #         add = "boxplot", add.params = list(fill = "white"))+
  #  stat_compare_means(label = "p.signif")+ # Add significance levels
  #  stat_compare_means(label.y = 50)
  #median_time <- fulltimeevo %>%
  
  
  listofsimulations <-
    sample(1:max(fulltimeevo$simulation_number), subset_n)
  fulltimeevo <-
    fulltimeevo %>% filter(simulation_number %in% listofsimulations)
  
  fulltimeevo$event <-
    factor(fulltimeevo$event,
           levels = c("Ongoing", "Expansion", "Extinction"))
  
  p002 <-
    ggplot(data = fulltimeevo, aes(y = time, x = event, fill = event)) + geom_violin(aes(alpha =
                                                                                           0.5)) +
    scale_alpha(guide = 'none')  +
    theme_classic2() + ylim(0, maxt) + scale_fill_brewer(palette = "Dark2", drop =
                                                           F) +
    stat_summary(
      fun.y = mean,
      geom = "point",
      shape = 23,
      size = 2
    ) +
    stat_summary(
      fun.y = median,
      geom = "point",
      shape = 20,
      size = 2
    ) +
    stat_summary(
      aes(label = round(..y.., 2)),
      fun.y = median,
      geom = "text",
      hjust = -1
    )
  
  p004 <-
    ggplot(data = fulltimeevo, aes(x = event, fill = event)) + geom_bar(aes(alpha =
                                                                              0.5)) +
    theme_classic2() + ylim(0, subset_n) + scale_fill_brewer(palette = "Dark2", drop =
                                                               F) +
    geom_text(aes(label = ..count..),
              stat = "count",
              position = position_stack())
  
  p003 <-
    ggplot(data = fulltimeevo, aes(fill = event, x = time))  + theme_classic2() +
    geom_dotplot(
      binaxis = 'x',
      stackdir = 'up',
      position = "dodge",
      method = "histodot",
      dotsize = .3,
      aes(fill = event)
    ) + scale_fill_brewer(palette = "Dark2", drop = F) + xlim(0, maxt)
  
  
  p005 <-
    ggarrange(
      ggarrange(
        p002,
        p004,
        nrow = 1,
        labels = c("A", "B"),
        legend = "none"
      ),
      p003,
      nrow = 2,
      labels = c("", "C"),
      common.legend = T
    )
  
  annotate_figure(p005, top = text_grob(
    paste("Fig", maintitle, sep = " "),
    color = "black",
    face = "bold",
    size = 14
  ))
  
}

##Plot global, driver and immune dnds for different simulations
plot_threeV4 = function (x,
                         subset_n = 20,
                         maxt = 100,
                         maintitle = "Summary simulation")
{
  x$simulation <- as.integer(x$simulation)
  listofsimulations <- sample(1:max(x$simulation), subset_n)
  x <- x %>% filter(simulation %in% listofsimulations)
  
  x$simulation <- as.factor(x$simulation)
  nb.cols <- length(unique(x$simulation))
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
  names(mycolors) = unique(x$simulation)
  
  #### Population count
  p1 <-
    ggplot(x, aes(
      x = t_0,
      y = count,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    labs(title = "Population size", y = 'Pop size', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") +  geom_abline(slope = 0,
                                                                     intercept = log10(16),
                                                                     col = "red") +
    scale_color_manual(values = mycolors) + xlim(0, maxt) + scale_y_log10()
  
  #### Mutation count
  p2 <-
    ggplot(x, aes(
      x = t_0,
      y = sum_na,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    labs(title = "Mutations", y = 'Non-silent muts', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") + scale_y_log10() + xlim(0, maxt) +
    scale_color_manual(values = mycolors)
  
  
  ###Global dN/dS
  p3 <-
    ggplot(x, aes(
      x = t_0,
      y = dnds,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    geom_abline(slope = 0,
                intercept = 0,
                col = "red") +
    labs(title = "Overall selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") + scale_y_log10() + xlim(0, maxt) +
    scale_color_manual(values = mycolors)
  
  ####Driver dN/dS
  dfx <-
    x %>% mutate(
      dnds_d = ifelse(is.na(dnds_d), 0, dnds_d),
      dnds_d = ifelse(is.infinite(dnds_d), count / 16, dnds_d)
    )
  if (nrow(dfx) > 0) {
    p4 <-
      ggplot(dfx, aes(
        x = t_0,
        y = dnds_d,
        color = simulation,
        alpha = 0.5
      )) +  geom_point() +
      geom_line() +
      geom_abline(slope = 0,
                  intercept = 0,
                  col = "red") +
      labs(title = "Driver selection", y = 'dN/dS Driver', x = "Generation") +
      theme_minimal() + theme(legend.position = "none") + scale_y_log10() +
      scale_color_manual(values = mycolors) + xlim(0, maxt)
  } else{
    p4 <- ggplot()
  }
  
  ####Immune dN/dS
  dfi <- x %>% filter(!is.na(dnds_i))
  if (nrow(dfi) > 0) {
    p5 <-
      ggplot(
        dfi %>% filter(!is.na(dnds_i), !is.infinite(dnds_i)),
        aes(
          x = t_0,
          y = dnds_i,
          color = simulation,
          alpha = 0.5
        )
      ) +  geom_point() +
      geom_line() +
      geom_abline(slope = 0,
                  intercept = 0,
                  col = "red") +
      labs(title = "Immune selection", y = 'dN/dS Immuno', x = "Generation") +
      theme_minimal() + theme(legend.position = "none") + scale_y_log10() +
      scale_color_manual(values = mycolors) + xlim(0, maxt)
  } else{
    p5 <- ggplot()
  }
  
  p6 <-
    ggplot(x, aes(
      x = t_0,
      y = numclones,
      color = simulation,
      alpha = 0.5
    )) +  geom_point() +
    geom_line() +
    labs(title = "Clone heterogeneity", y = 'Number of clones', x = "Generation") +
    theme_minimal() + theme(legend.position = "none") +
    scale_color_manual(values = mycolors) + xlim(0, maxt)
  
  finalp1 <-
    ggarrange(p1,
              p2,
              p6,
              p3,
              p4,
              p5,
              labels = "AUTO",
              common.legend = F)
  
  annotate_figure(finalp1, top = text_grob(
    paste("Fig", maintitle, sep = " "),
    color = "black",
    face = "bold",
    size = 14
  ))
}

##Plot metrics from different classes
plot_metric_class = function (x,
                              subset_n = 80,
                              maxt = 100,
                              maintitle = "Summary simulation")
{
  #library(viridis)
  library("ggsci")
  library("ggplot2")
  library("gridExtra")
  
  listofsimulations <-
    sample(1:max(x$simulation), subset_n, replace = T)
  x <- x %>% filter(simulation %in% listofsimulations)
  
  
  medians <- aggregate(dnds ~  t_0 + class, x, median)
  #medians <- round(medians)
  
  x$simulation <- as.factor(x$simulation)
  x$class <- as.factor(x$class)
  
  #Colors
  #nb.cols <- length(unique(x$class))
  #mycolors <- colorRampPalette(brewer.pal(6, "Purples"))(nb.cols)
  #names(mycolors) = unique(x$class)
  
  #mycolors <- c("MSS" = "grey", "MSI" = "blue", "POLE" = "black")
  
  p1 <-
    ggplot(x  %>% filter(t_0 <= maxt), aes(
      x = factor(t_0),
      y = dnds,
      fill = factor(class)
    )) +  #geom_point(alpha = 0.5, size = 1) +
    geom_boxplot(position = "dodge2",
                 outlier.shape = NA,
                 alpha = 0.5) +
    #geom_line() +
    geom_abline(slope = 0,
                intercept = 0,
                col = "red") +
    labs(title = "Overall selection over time", y = 'dN/dS', x = "Generation") +
    theme_minimal()  + scale_y_log10()  +
    #geom_smooth(method = "loess", alpha = 0.05, size = 1, span = 1) +
    scale_fill_npg() #+
  #scale_color_manual(values = mycolors, drop =F)
  
  p2 <-
    ggplot(x %>% filter(t_0 <= maxt),
           aes(
             x = factor(t_0),
             y = count,
             fill = factor(class)
           )) +  #geom_point(alpha = 0.5, size = 1) +
    geom_boxplot(position = "dodge2",
                 outlier.shape = NA,
                 alpha = 0.5) + #geom_line() +
    labs(title = "Population size", y = 'Pop size', x = "Generation") +
    theme_minimal()  + scale_y_log10() +
    geom_abline(slope = 0,
                intercept = log10(16),
                col = "red") +
    scale_fill_npg() #+
  # xlim(0,maxt)
  #scale_color_manual(values = mycolors, drop =F)
  
  #ggplot(fulltime_neutral,aes(x = factor(t_0),y = count, fill=factor(class))) +
  #  geom_boxplot(position="dodge2",outlier.shape = NA)
  
  p3 <-
    ggplot(x %>% filter(t_0 <= maxt),
           aes(
             x = factor(t_0),
             y = sum_na,
             fill = factor(class)
           )) +  #geom_point(alpha = 0.5, size = 1) +
    geom_boxplot(position = "dodge2",
                 outlier.shape = NA,
                 alpha = 0.5) +
    #geom_line() +
    labs(title = "Mutations", y = 'Non-silent muts', x = "Generation") +
    theme_minimal() + scale_y_log10() +
    scale_fill_npg()
  #scale_color_manual(values = mycolors, drop =F)
  
  finalp1 <-
    ggarrange(
      p2,
      p3,
      p1,
      nrow = 3,
      labels = "AUTO",
      common.legend = T,
      legend = "bottom"
    )
  
  annotate_figure(finalp1, top = text_grob(
    paste("Fig", maintitle, sep = " "),
    color = "black",
    face = "bold",
    size = 14
  ))
  
}

plot_metric_classevo = function (fulltimeevo,
                                 maxt = 100,
                                 maintitle = "Summary simulation")
{
  library(ggsci)
  #listofsimulations<- sample(1:max(fulltimeevo$simulation_number),subset_n)
  #fulltimeevo<-fulltimeevo %>% filter(simulation_number %in% listofsimulations)
  
  fulltimeevo$event <-
    factor(fulltimeevo$event,
           levels = c("Ongoing", "Expansion", "Extinction"))
  fulltimeevo$class <- as.factor(fulltimeevo$class)
  
  p002 <-
    ggplot(data = fulltimeevo %>% filter(event != "Ongoing"),
           aes(y = time, x = class, fill = class)) + geom_violin(aes(alpha = 0.5)) +
    scale_alpha(guide = 'none')  +
    #stat_summary(fun.y=mean, geom="point", shape=23, size=2, position=position_nudge(x = 0.1)) +
    stat_summary(
      fun.y = median,
      geom = "point",
      shape = 20,
      size = 2,
      position = "dodge2"
    ) +
    stat_summary(
      aes(label = round(..y.., 2)),
      fun.y = median,
      geom = "text",
      hjust = -1,
      position = "dodge2"
    ) +
    theme_classic2() + ylim(0, maxt) + scale_color_npg(alpha = 0.5) + facet_wrap( ~
                                                                                    event) + rotate_x_text(angle = 45)
  
  p004 <-
    ggplot(data = fulltimeevo, aes(x = class, fill = class)) + geom_bar(
      aes(
        label = ..count..,
        alpha = 0.5,
        fill = event
      ),
      stat = "count",
      position = position_fill()
    ) +
    theme_classic2() + scale_fill_npg() + rotate_x_text(angle = 45)#+ scale_color_nejm()
  #geom_text(aes(label=..count..),stat="count",position=position_stack()) + facet_wrap(~event)
  
  p003 <-
    ggplot(data = fulltimeevo, aes(fill = event, x = time, alpha = 0.5))  + theme_classic2() +
    geom_dotplot(
      binaxis = 'x',
      stackdir = 'up',
      position = "dodge",
      method = "histodot",
      #stackgroups = T,
      dotsize = .9,
      aes(fill = event)
    ) + scale_fill_npg() + facet_wrap( ~ class) + guides(alpha = F) + rotate_x_text(angle = 45)
  
  
  p005 <-
    ggarrange(
      ggarrange(
        p004,
        p002,
        nrow = 1,
        labels = c("A", "B"),
        legend = "none"
      ),
      p003,
      nrow = 2,
      labels = c("", "C"),
      common.legend = T
    )
  
  annotate_figure(p005, top = text_grob(
    paste("Fig", maintitle, sep = " "),
    color = "black",
    face = "bold",
    size = 14
  ))
  #p004
}

plot_metric_classevo2 = function (fulltimeevo,
                                  maxt = 100,
                                  maintitle = "Summary simulation")
{
  library(ggsci)
  
  fulltimeevo$event <-
    factor(fulltimeevo$event,
           levels = c("Ongoing", "Expansion", "Extinction"))
  
  p004 <- ggboxplot(
    data = fulltimeevo %>% filter(event != "Ongoing"),
    y = "time",
    x = "class",
    # Add regression line
    fill = "event"#,                                  # Add confidence interval
    #add="jitter"
  ) + theme_classic2() + ylim(0, maxt) + scale_fill_npg(alpha = 0.5)  +
    facet_wrap( ~ event)  + rotate_x_text(angle = 45) #+ stat_compare_means()
  #+ #geom_violin(aes(alpha=0.5)) +
  #scale_alpha(guide = 'none')
  #stat_summary(fun.y=mean, geom="point", shape=23, size=2, position=position_nudge(x = 0.1)) +
  #stat_summary(fun.y=median, geom="point", shape=20, size=2, position="dodge2") +
  #stat_summary(aes(label=round(..y..,2)),fun.y=median, geom="text",hjust = -1, position="dodge2") +
  #theme_classic2() + ylim(0,maxt) + scale_color_npg(alpha = 0.5) + facet_wrap(~event) + rotate_x_text(angle= 45)
  
  
  p002 <-
    ggplot(data = fulltimeevo %>% filter(event != "Ongoing"),
           aes(x = as.factor(class), fill = class)) +
    geom_bar(
      aes(
        label = ..count..,
        alpha = 0.5,
        fill = event
      ),
      stat = "count",
      position = position_fill()
    ) +
    theme_classic2() + scale_fill_npg() + rotate_x_text(angle = 45)#+ scale_color_nejm()
  #geom_text(aes(label=..count..),stat="count",position=position_stack()) + facet_wrap(~event)
  
  #p003<-ggplot(data=fulltimeevo, aes(fill=event, x = time, alpha =0.5))  + theme_classic2() +
  #  geom_dotplot(binaxis='x',
  #               stackdir = 'up',
  #               position = "dodge",
  #               method = "histodot",
  #               #stackgroups = T,
  #               dotsize = .9,
  #               aes(fill=event)) + scale_fill_nejm() + facet_wrap(~class) +guides(alpha=F) + rotate_x_text(angle= 45)
  
  
  p005 <-
    ggarrange(
      p002,
      p004,
      nrow = 1,
      labels = c("A", "B"),
      legend = "none",
      common.legend = T
    )
  
  annotate_figure(p005, top = text_grob(
    paste("Fig", maintitle, sep = " "),
    color = "black",
    face = "bold",
    size = 14
  ))
  #p004
}


plot_statistic_clasevo1 = function(x)
{
  library()
  p1 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(event != "Ongoing", event != "Extinction"),
      x = class2,
      y = time,
      #grouping.var =  event,
      type = "np",
      #plot.type = "violin",
      pairwise.comparisons = T,
      bf.message = F,
      results.subtitle = F,
      messages = F,
      mean.size = 1
      
    ) + rotate_x_text(45)
  
  
  p2 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(fin_dnds != Inf, event != "Ongoing", event != "Extinction"),
      x = class2,
      y = fin_dnds,
      #grouping.var =  event,
      type = "np",
      #plot.type = "violin",
      pairwise.comparisons = T,
      bf.message = F,
      results.subtitle = F,
      messages = F,
      mean.size = 1
    ) + rotate_x_text(45) + scale_y_log10()
  
  p3 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(fin_dndsd != Inf, event != "Ongoing", event != "Extinction"),
      x = class2,
      y = fin_dndsd,
      #grouping.var =  event,
      type = "np",
      #plot.type = "violin",
      pairwise.comparisons = T,
      bf.message = F,
      results.subtitle = F,
      messages = F,
      mean.size = 1
    ) + rotate_x_text(45) + scale_y_log10()
  
  p4 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(fin_dndsi != Inf, event != "Ongoing", event != "Extinction"),
      x = class2,
      y = fin_dndsi,
      #grouping.var =  event,
      type = "np",
      plot.type = "violin",
      pairwise.comparisons = T,
      bf.message = F,
      results.subtitle = F,
      messages = F,
      mean.size = 1
    ) + rotate_x_text(45) + scale_y_log10()
  
  
  ggarrange(p1, p2, p3, p4, labels = "AUTO")
  
}


plot_statistic_clasevo2 = function(x)
{
  p2 <- ggstatsplot::ggbarstats(
    data = x ,
    main = event,
    condition = class2,
    sampling.plan = "jointMulti",
    perc.k = 1,
    x.axis.orientation = "slant",
    #ggtheme = hrbrthemes::theme_modern_rc(),
    ggstatsplot.layer = FALSE,
    ggplot.component = ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic")),
    palette = "nrc_npg",
    package = "ggsci",
    messages = FALSE,
    bf.message = F,
    results.subtitle = F
  ) + rotate_x_text(45)
  
  p3 <- ggstatsplot::ggbarstats(
    data = x %>% mutate(class2 = ifelse(
      class == "I:0.99",
      "HIGH",
      ifelse(
        class == "I:0.95",
        "HIGH",
        ifelse(
          class == "I:0.75",
          "HIGH",
          ifelse(class == "I:0.50", "MID", "LOW")
        )
      )
    )),
    main = class2,
    condition = event,
    sampling.plan = "jointMulti",
    perc.k = 1,
    x.axis.orientation = "slant",
    #ggtheme = hrbrthemes::theme_modern_rc(),
    ggstatsplot.layer = FALSE,
    ggplot.component = ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic")),
    palette = "default_npg",
    package = "ggsci",
    messages = FALSE,
    bf.message = F,
    results.subtitle = F
  ) + rotate_x_text(45)
  
  ggarrange(p2, p3, labels = "AUTO")
  
}

plot_statistic_clasevo3 = function(x)
{
  p01 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(event == "Expansion"),
      x = class,
      y = time,
      plot.type = "violin",
      type = "np",
      bf.message = F,
      results.subtitle = F,
      #notch = T,
      linetype = "dashed",
      outlier.shape = NA,
      messages = F,
      p.adjust.method = "BH",
      package = "ggsci",
      palette = "nrc_npg",
      mean.plotting = T,
      mean.label.size = 2,
      mean.color = "darkgrey",
      mean.size = 1,
      sample.size.label = F,
      pairwise.comparisons = T
    ) + rotate_x_text(45)
  
  
  p02 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(event == "Extinction"),
      x = class,
      y = time,
      plot.type = "violin",
      type = "np",
      bf.message = F,
      results.subtitle = F,
      #notch = T,
      linetype = "dashed",
      outlier.shape = NA,
      messages = F,
      p.adjust.method = "BH",
      package = "ggsci",
      palette = "nrc_npg",
      mean.plotting = T,
      mean.label.size = 2,
      mean.color = "darkgrey",
      mean.size = 1,
      sample.size.label = F,
      pairwise.comparisons = T
    ) + rotate_x_text(45)
  
  ggarrange(p01, p02, labels = "AUTO")
  
}

plot_statistic_clasevo4 = function(x)
{
  p03 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(
        event != "Ongoing",
        !is.infinite(fin_dnds),
        event != "Extinction"
      ),
      x = class,
      y = fin_dnds,
      plot.type = "box",
      type = "np",
      bf.message = F,
      results.subtitle = F,
      pairwise.comparisons = T,
      pairwise.annotation = "asterisk",
      pairwise.display = "all",
      notch = T,
      #linetype = "dashed",
      outlier.shape = NA,
      messages = T,
      p.adjust.method = "BH",
      package = "ggsci",
      palette = "nrc_npg",
      ggtheme = theme_classic(),
      mean.plotting = F,
      mean.label.size = 1,
      mean.color = "darkgrey",
      ylab = "Global dN/dS",
      xlab = "P(Immune Attack)",
      mean.size = 1,
      sample.size.label = F
    ) + rotate_x_text(45) + geom_hline(yintercept = 1,
                                       linetype = "dashed",
                                       color = "red")
  
  
  p04 <-
    ggstatsplot::ggbetweenstats(
      data = x %>% filter(
        event != "Ongoing",
        !is.infinite(fin_dndsd),
        event != "Extinction"
      ),
      x = class,
      y =  fin_dndsd,
      plot.type = "box",
      type = "np",
      bf.message = F,
      results.subtitle = F,
      pairwise.comparisons = T,
      pairwise.annotation = "asterisk",
      pairwise.display = "all",
      #notch = T,
      notch = F,
      #linetype = "dashed",
      ggtheme = theme_classic(),
      outlier.shape = NA,
      messages = F,
      p.adjust.method = "BH",
      package = "ggsci",
      palette = "nrc_npg",
      mean.plotting = F,
      mean.label.size = 1,
      mean.color = "darkgrey",
      ylab = "Driver dN/dS",
      xlab = "P(Immune Attack)",
      mean.size = 1,
      sample.size.label = F
    ) + rotate_x_text(45) + coord_trans(y = "log") + geom_hline(yintercept =
                                                                  1,
                                                                linetype = "dashed",
                                                                color = "red")
  
  
  
  p05 <-
    ggstatsplot::ggbetweenstats(
      data =  x %>% filter(
        event != "Ongoing",
        !is.infinite(fin_dndsi),
        event != "Extinction"
      ),
      x = class,
      y =  fin_dndsi,
      plot.type = "box",
      type = "np",
      bf.message = F,
      results.subtitle = F,
      pairwise.comparisons = T,
      pairwise.annotation = "asterisk",
      pairwise.display = "all",
      ggtheme = theme_classic(),
      #notch = T,
      notch = F,
      #linetype = "dashed",
      outlier.shape = NA,
      messages = F,
      p.adjust.method = "BH",
      package = "ggsci",
      palette = "nrc_npg",
      mean.plotting = F,
      mean.label.size = 1,
      mean.color = "darkgrey",
      ylab = "Immune dN/dS",
      xlab = "P(Immune Attack)",
      mean.size = 1,
      sample.size.label = F
    ) + rotate_x_text(45) + coord_trans(y = "log") + geom_hline(yintercept =
                                                                  1,
                                                                linetype = "dashed",
                                                                color = "red")
  
  ggarrange(p03, p04, p05, nrow = 1, labels = "AUTO")
}

#get number of clones at time t
get_numberclones = function(M, t) 
{
  length(M %>% filter(t_0 == !!t) %>% pull(clone_id) %>% unique())
}
# Return the clone tree with the ggmuller naming convention
get_clone_tree = function(M)
{
  M %>%
    dplyr::select(parent_clone_id, clone_id) %>%
    distinct() %>%
    rename(Parent = parent_clone_id, Identity = clone_id) %>%
    filter(Parent + Identity > 0) %>%
    mutate(Parent = paste(Parent), Identity = paste(Identity))
}

## get psurv as fitness value
get_clone_fitness = function(M)
{
  M %>%
    dplyr::select(psurv, clone_id) %>%
    group_by(clone_id) %>%
    distinct() %>%
    rename(Identity = clone_id,
           Fitness = psurv) %>%
    ungroup() %>%
    mutate(Identity = paste(Identity),
           Fitness = Fitness)
}

##Create muller plot and add clone tree
mullerplot = function(M,
                      t_from = 0,
                      t_to = max(M$t_0),
                      ps_0 = 0.5,
                      cutoff = 0,
                      legend = FALSE,
                      add_clone_tree = TRUE,
                      ...) {
  require(ggpubr)
  require(ggmuller)
  
  pio::pioTit("Muller plot for the simulation")
  
  nC = length(unique(M$clone_id))
  
  pio::pioStr("All simulation stats", '\n')
  pio::pioStr("Number of clones: ", paste0(nC, ' - initial clone with id 0\n'))
  if (nC > 20 &
      legend)
    message("More than 20 clones, if the legend is messy use `legend = FALSE`.")
  
  # Reduce time range
  pio::pioStr("Using time range for plot: ", paste0(t_from, ' - ', t_to, '\n'))
  # M = M %>% filter(t_0 >= t_from & t_0 <= t_to)
  
  # Clone tree
  Clone_phylogeny = get_clone_tree(M %>% filter(t_0 >= t_from &
                                                  t_0 <= t_to))
  # %>%
  #   filter(Parent != Identity)
  #
  # Population counts
  Pop_counts = M %>%
    filter(t_0 >= t_from & t_0 <= t_to) %>%
    rename(Generation = t_0, Identity = clone_id) %>%
    group_by(Generation, Identity) %>%
    summarise(Counts_Population = n()) %>%
    ungroup()
  
  # Spread the dataframe to fill in automatically all missing 0s
  Pop_counts = Pop_counts %>%
    spread(Identity, Counts_Population) %>%
    replace(., is.na(.), 0)
  
  Pop_counts = Pop_counts %>%
    reshape2::melt(id = "Generation") %>%
    rename(Identity = variable, Population =  value) %>%
    as_tibble() %>%
    mutate(Identity = paste(Identity))
  
  ###Fitness
  Fitness = get_clone_fitness(M %>% filter(t_0 >= t_from &
                                             t_0 <= t_to))
  
  Pop_counts = Pop_counts %>% left_join(Fitness, by = "Identity")
  
  Muller_df <-
    get_Muller_df(Clone_phylogeny, Pop_counts, cutoff = cutoff)
  
  leg = theme(legend.position = 'bottom',
              legend.key.size = unit(.3, "cm"))
  
  p1 = Muller_plot(Muller_df, ...) +
    labs(title = "Clone proportions over time") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")
  
  p2 = Muller_pop_plot(Muller_df, ...) +
    labs(title = "Clone size over time") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")
  
  p3 = Muller_plot(Muller_df, colour_by = "Fitness", ...) +
    scale_fill_gradientn(colors = c('steelblue', 'darkred'),
                         breaks = seq(0, 1, 0.1)) +
    #
    # scale_fill_gradientn(colors = c('steelblue', 'gray', 'darkred'), breaks = c(0, ps_0, 1)) +
    labs(title = "Fitness: psurv", y = "Clone proportions") +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")
  
  p5 = dnds_plot(M, ...) +
    labs(title = "Overall Selection", y = "dN/dS")
  
  if (legend)
  {
    p1 = p1 + leg
    p2 = p2 + leg
    p3 = p3 + leg
    # p4 = p4 + leg
    p5 = p5 + leg
  }
  
  figure = NULL
  
  # Add clone tree
  cl_plot = NULL
  if (add_clone_tree)
  {
    cl_plot = clonetree_plot(M, ps_0 = ps_0 , t_to = t_to)
    figure = ggarrange(p2, p1, p3, p5, cl_plot, nrow = 1, ncol = 5)
  }
  else
    figure = ggarrange(p2, p1, p3, p5, nrow = 1, ncol = 4)
  figure
}

##Create clone tree plot
clonetree_plot = function(M,
                          ps_0 = 0.5,
                          t_to = max(M$t_0),
                          ...) {
  require(tidygraph)
  require(ggraph)
  
  pio::pioTit("Clone tree plot for the simulation")
  
  alive_clones = M %>%
    filter(alive) %>%
    pull(clone_id) %>%
    unique()
  
  # Reduce time range
  pio::pioStr("Using time range: t <= ", t_to, '\n')
  M = M %>% filter(t_0 <= t_to)
  
  # Clone tree
  Clone_phylogeny = get_clone_tree(M)
  
  print(Clone_phylogeny)
  
  # tree <- adj_matrix_to_tree(Clone_phylogeny)
  #
  # library(ape)
  # tree$tip.label <- 1:length(tree$tip.label) # optional
  # tree$node.label <- (length(tree$tip.label) + 1):10 # optional
  #
  # plot(tree, show.node.label = TRUE, show.tip.label = TRUE, tip.color = "red")
  # lC = unique(unlist(Clone_phylogeny))
  # nC = length(unique(lC))
  
  # Use a tidygraph
  tb_tree = as_tbl_graph(Clone_phylogeny)
  
  # Attributes
  attributes = get_clone_fitness(M)
  attributes$alive = FALSE
  
  attributes = attributes %>%
    mutate(alive = ifelse(Identity %in% alive_clones, TRUE, FALSE))
  
  tb_tree = tb_tree %>%
    activate(nodes) %>%
    mutate(Identity = name) %>%
    left_join(attributes, by = 'Identity')
  
  # Plot call
  ggraph(tb_tree, 'tree') +
    geom_edge_diagonal(
      arrow = arrow(length = unit(2, 'mm')),
      end_cap = circle(5, 'mm'),
      start_cap  = circle(5, 'mm')
    ) +
    geom_node_point(aes(colour = Fitness, shape = alive), size = 8) +
    geom_node_text(aes(label = Identity),
                   colour = 'black', vjust = 0.4) +
    theme_void(base_size = 8) +
    theme(legend.position = 'bottom',
          legend.key.size = unit(.3, "cm")) +
    scale_color_gradientn(colors = c('steelblue', 'darkred'),
                          breaks = seq(0, 1, 0.1)) +
    # scale_color_gradientn(colors = c('steelblue', 'gray', 'darkred'), breaks = c(0, ps_0, 1)) +
    #
    # scale_color_distiller(palette = 'RdBu', direction = -1) +
    labs(title = paste("Clone tree for t <=", t_to),
         subtitle = "Clones coloured by fitness, alive clone status annotated for the end of the simulation")
}

##Create basic dnds plot
dnds_plot = function(M, t_to = max(M$t_0), ...) {
  M = M %>% filter(t_0 <= t_to)
  M = M %>% filter(strategy != "Death")
  ####Get the population parameters over time for the population size
  ####Get the population parameters over time for the population size
  M = M %>% mutate(t_0 = ifelse(is.na(t_0), max(M$t_0) + 1, t_0))
  
  pop = M %>%
    group_by(t_0) %>%
    summarise(count = n()) %>%
    mutate(expon = 2 ^ (t_0))
  
  ####Get the population parameters for the population dN/dS
  DNDS = M %>%
    group_by(t_0) %>%
    summarise(
      gnad = sum(gnad),
      gnap = sum(gnap),
      gnai = sum(gnai),
      gnae = sum(gnae),
      gnak = sum(gnak),
      gns = sum(gns),
      gnsi = sum(gnsi),
      gnsd = sum(gnsd),
      median_dnds = median(dnds),
      median_dndsi = median(dnds_i),
      median_dndsd = median(dnds_d)
    ) %>%
    mutate(
      dnds = compute_gdnds(gnad, gnap, gnai, gnae, gnak, gns, gnsi, gnsd),
      dnds_i = compute_gdnds_i(gnai, gnsi),
      dnds_d = compute_gdnds_d(gnad, gnsd),
      dnds_lci =  dnds - (1.96 * (sqrt(
        dnds / (gnad + gnap + gnai + gnae + gnak + gns + gnsi)
      ))),
      dnds_hci =  dnds + (1.96 * (sqrt(
        dnds / (gnad + gnap + gnai + gnae + gnak + gns + gnsi)
      ))),
      dnds_ilci = dnds_i - (1.96 * (sqrt(
        dnds_i / (gnai + gnsi)
      ))),
      dnds_ihci = dnds_i + (1.96 * (sqrt(
        dnds_i / (gnai + gnsi)
      ))),
      dnds_dlci = dnds_d - (1.96 * (sqrt(
        dnds_d / (gnad + gnsd)
      ))),
      dnds_dhci = dnds_d + (1.96 * (sqrt(
        dnds_d / (gnad + gnsd)
      )))
    )
  
  ##Get the totals
  totals = M %>%
    group_by(t_0) %>%
    summarise(
      gnad = sum(gnad),
      gnap = sum(gnap),
      gnai = sum(gnai),
      gnae = sum(gnae),
      gnak = sum(gnak),
      gns = sum(gns),
      gnsi = sum(gnsi),
      gnsd = sum(gnsd),
      sum_na = gnad + gnap + gnai + gnae + gnak,
      sum_ns = gns + gnsi + gnad
    )
  
  ###In order to plot we need to join a simplified table
  time_series = full_join(pop, DNDS) %>% full_join(totals) %>%
    dplyr::select(
      t_0,
      median_dnds,
      median_dndsi,
      median_dndsd,
      dnds,
      dnds_lci,
      dnds_hci,
      dnds_i,
      dnds_ilci,
      dnds_ihci,
      dnds_d,
      dnds_dlci,
      dnds_dhci,
      count,
      expon,
      sum_na,
      sum_ns,
      gnad,
      gnap,
      gnai,
      gnae,
      gnak,
      gns,
      gnsi,
      gnsd
    )
  
  #time_series = time_series %>% filter(t_0 <= 26)
  
  
  #quartz()
  p1 <- ggplot(
    time_series %>%
      filter(!is.na(dnds), dnds != Inf),
    aes(
      x = t_0,
      y = dnds,
      ymax = dnds_hci,
      ymin = dnds_lci
    )
  ) +
    geom_point() +
    geom_pointrange(col = "plum4") +
    geom_line(aes(col = "#c5b2c5", alpha = 0.2)) +
    cowplot::theme_cowplot() +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") +
    #ylim(c(0,max(time_series$dnds))) +
    xlim(c(0, max(time_series$t_0))) +
    labs(title = "Overall Selection", y = 'dN/dS', x = "Generation")
  
  p2 <- ggplot(
    time_series %>%
      filter(!is.na(dnds_i), dnds_i != Inf),
    aes(
      x = t_0 ,
      y = dnds_i,
      ymax = dnds_ihci,
      ymin = dnds_ilci
    )
  ) +
    geom_point() +
    geom_pointrange(col = "plum4") +
    geom_line(aes(col = "#c5b2c5", alpha = 0.2)) +
    cowplot::theme_cowplot() +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") +
    #ylim(c(0,max(time_series$dnds_i))) +
    xlim(c(0, max(time_series$t_0))) +
    labs(title = "Immune Selection", y = 'dN/dS immunopeptidome', x = "Generation")
  
  p3 <- ggplot(
    time_series %>%
      filter(!is.na(dnds_d), dnds_d != Inf),
    aes(
      x = t_0 ,
      y = dnds_d,
      ymax = dnds_dhci,
      ymin = dnds_dlci
    )
  ) +
    geom_point() +
    geom_pointrange(col = "plum4") +
    geom_line(aes(col = "#c5b2c5", alpha = 0.2)) +
    cowplot::theme_cowplot() +
    geom_abline(slope = 0,
                intercept = 1,
                col = "red") +
    #ylim(c(0,max(time_series$dnds_i))) +
    xlim(c(0, max(time_series$t_0))) +
    labs(title = "Driver Selection", y = 'dN/dS drivers', x = "Generation")
  
  if (nrow(time_series %>% filter(!is.na(dnds_i), dnds_i != Inf)) == 0 &
      nrow(time_series %>% filter(!is.na(dnds_d), dnds_d != Inf)) == 0) {
    ggarrange(p1, legend = "none")
  }
  else if (nrow(time_series %>% filter(!is.na(dnds_i), dnds_i != Inf)) == 0 &
           nrow(time_series %>% filter(!is.na(dnds_d), dnds_d != Inf)) != 0) {
    ggarrange(p1, p3, nrow = 2, legend = "none")
  }
  else if (nrow(time_series %>% filter(!is.na(dnds_i), dnds_i != Inf)) != 0 &
           nrow(time_series %>% filter(!is.na(dnds_d), dnds_d != Inf)) == 0) {
    ggarrange(p1, p2, nrow = 2, legend = "none")
  }
  else {
    ggarrange(p1, p3, p2, nrow = 3, legend = "none")
  }
}

##Create basic dnds plot per clone id
plot_clone_dnds = function(x) {
  library(RColorBrewer)
  
  data = x %>% filter(strategy != "Death") %>%
    group_by(clone_id, t_0) %>%
    summarise(
      gnad = sum(gnad),
      gnap = sum(gnap),
      gnai = sum(gnai),
      gnae = sum(gnae),
      gnak = sum(gnak),
      gns = sum(gns),
      gnsi = sum(gnsi),
      gnsd = sum(gnsd),
      count = n()
    ) %>%
    mutate(
      dnds = compute_gdnds(gnad, gnap, gnai, gnae, gnak, gns, gnsi, gnsd),
      dnds_i = compute_gdnds_i(gnai, gnsi),
      dnds_d = compute_gdnds_d(gnad, gnsd),
      escape = ifelse(
        (gnae / count) > 0.99 ,
        "clonal",
        ifelse(gnae == 0, "absent", "subclonal")
      )
    ) %>%
    ungroup()
  
  life_clone_length = x %>%
    group_by(clone_id) %>%
    summarise(L = max(t_0, na.rm = T) - min(t_0, na.rm = T)) %>%
    ungroup() %>%
    # arrange(desc(L))
    filter(L > 10) %>% pull(clone_id)
  
  clone_size = x %>%
    group_by(clone_id, t_0) %>%
    summarise(cells = n()) %>%
    ungroup() %>%
    group_by(t_0) %>%
    mutate(N = sum(cells),
           proportion = cells / N) %>%
    ungroup()
  
  data = data %>%
    filter(clone_id %in% life_clone_length)
  
  clone_size = clone_size  %>%
    filter(clone_id %in% life_clone_length)
  
  data = data %>%
    left_join(clone_size, by = c('clone_id', 't_0'))
  
  #Fix oreder of levels for escape
  data$escape = factor(data$escape, levels = c("absent", "subclonal", "clonal"))
  
  #Fix dN/dS when values are infinite for Drivers only (We assume that the dN/dS cant be larger than the size of the population)
  data <-
    data %>% mutate(
      dnds_d = ifelse(is.infinite(dnds_d), cells, dnds_d),
      dnds_d = ifelse(is.na(dnds_d), 0, dnds_d)
    )
  
  #Colors
  nb.cols <- length(unique(data$clone_id))
  mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
  names(mycolors) = unique(data$clone_id)
  
  
  p1 <-  ggplot(data %>% filter(!is.na(dnds), dnds != Inf),
                aes(
                  x = t_0,
                  y = dnds,
                  color = paste(clone_id)
                )) +
    geom_line() +
    xlim(c(0, max(data$t_0))) +
    geom_point(aes(
      size = proportion,
      alpha = 0.8,
      shape = escape
    )) +
    cowplot::theme_cowplot() + scale_alpha(guide = 'none') + geom_abline(slope = 0,
                                                                         intercept = 0,
                                                                         col = "red")   +
    labs(title = "Overall Selection", y = 'dN/dS', x = "Generation") + scale_y_log10() +
    scale_shape_manual(
      name = "Escape",
      labels = c("absent", "subclonal", "clonal"),
      values = c(16, 17, 15),
      drop = F
    ) +
    scale_color_manual(values = mycolors, drop = F)
  
  p2 <-  ggplot(data %>% filter(!is.na(dnds_i), dnds_i != Inf),
                aes(
                  x = t_0,
                  y = dnds_i,
                  color = paste(clone_id)
                )) +
    geom_line() +
    xlim(c(0, max(data$t_0))) +
    geom_point(aes(
      size = proportion,
      alpha = 0.8,
      shape = escape
    )) +
    cowplot::theme_cowplot() + scale_alpha(guide = 'none')  + geom_abline(slope = 0,
                                                                          intercept = 0,
                                                                          col = "red")   +
    labs(title = "Immunopeptidome Selection", y = 'dN/dS', x = "Generation") + scale_y_log10() +
    scale_shape_manual(
      name = "Escape",
      labels = c("absent", "subclonal", "clonal"),
      values = c(16, 17, 15),
      drop = F
    ) +
    scale_color_manual(values = mycolors)
  
  
  
  p3 <-  ggplot(data %>% filter(!is.na(dnds_d), dnds_d != Inf),
                aes(
                  x = t_0,
                  y = dnds_d,
                  color = paste(clone_id)
                )) +
    geom_line() +
    xlim(c(0, max(data$t_0))) +
    geom_point(aes(
      size = proportion,
      alpha = 0.8,
      shape = escape
    )) +
    cowplot::theme_cowplot() + scale_alpha(guide = 'none') + geom_abline(slope = 0,
                                                                         intercept = 0,
                                                                         col = "red")   +
    labs(title = "Driver Selection", y = 'dN/dS', x = "Generation") + scale_y_log10() +
    scale_shape_manual(
      name = "Escape",
      labels = c("absent", "subclonal", "clonal"),
      values = c(16, 17, 15),
      drop = F
    ) +
    scale_color_manual(values = mycolors)
  
  p4 <- ggplot(data,
               aes(
                 x = t_0,
                 y = count,
                 color = paste(clone_id)
               )) +
    geom_line() +
    xlim(c(0, max(data$t_0))) +
    geom_point(aes(
      size = proportion,
      alpha = 0.8,
      shape = escape
    )) +
    cowplot::theme_cowplot() + scale_alpha(guide = 'none') + geom_abline(slope = 0,
                                                                         intercept = log10(16),
                                                                         col = "red")   +
    labs(title = "Clonal populations", y = 'Size', x = "Generation") + scale_y_log10() +
    scale_shape_manual(
      name = "Escape",
      labels = c("absent", "subclonal", "clonal"),
      values = c(16, 17, 15),
      drop = F
    ) +
    scale_color_manual(values = mycolors)
  
  ggarrange(
    p1,
    p2,
    p3,
    p4,
    nrow = 4,
    labels = "AUTO",
    common.legend = T,
    legend = "bottom"
  )
  
}