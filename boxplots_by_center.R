
library(feather)
library(tidyverse)
library(Rtsne)
library(ggplot2)
library(plotly)
library(tidyverse)
library(htmlwidgets)
library(ggpubr)
library(rstatix)
library(dplyr)

setwd(dir = "D:/escritorio/tSNE_entropia_R")


  
  #### READ CSVs WITH SPECTRAL FEATURES ####
  n_Cal <- read.table(file = "nonorm_spectral_california_5s_no_overlapping.csv", sep = ",", header = TRUE)
  n_Med <- read.table(file = "nonorm_spectral_antioquia_5s_no_overlapping.csv", sep = ",", header = TRUE)
  n_Fin <- read.table(file = "nonorm_spectral_finland_5s_no_overlapping.csv", sep = ",", header = TRUE)
  n_Iow <- read.table(file = "nonorm_spectral_iowa_5s_no_overlapping.csv", sep = ",", header = TRUE)
  
  
  
  
  #MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
  Med_Par <- read.table(file = "participants_med.tsv", sep = "\t", header = TRUE)
  Med_Par <- Med_Par %>%
    separate(participant_id, c("Sub", "subject"), "-")
  Med_Par$subject <- as.numeric(Med_Par$subject)
  
  Med_Par <- Med_Par %>%
    dplyr::select(subject, group, age, gender)
  #  dplyr::select(subject, group, subgroup, moca_total, age)
  
  n_Med <- n_Med %>%
    inner_join(Med_Par, by="subject")
  n_Med$subject <- paste0("m-",n_Med$subject)
  
  
  
  #California Metadata
  Cal_par <- read.table(file = "participants_cal.tsv", sep = "\t", header = TRUE)
  Cal_par <- Cal_par %>%
    separate(participant_id, c("Sub", "subject"), "-")
  
  Cal_par <- Cal_par %>%
    dplyr::select(subject, age, gender)
  #  dplyr::select(subject, group, subgroup, moca_total, age)
  
  n_Cal <- n_Cal %>%
    inner_join(Cal_par, by="subject")
  
  n_Cal$subject <- paste0("c-",n_Cal$subject)
  
  
  
  
  #Finland metadata
  Fin_Par <- read.table(file = "participants_fin.tsv", sep = "\t", header = TRUE)
  Fin_Par <- Fin_Par %>%
    separate(participant_id, c("Sub", "subject"), "-")
  Fin_Par$subject <- as.numeric(Fin_Par$subject)
  
  Fin_Par <- Fin_Par %>%
    dplyr::select(subject, group, age, gender)
  #  dplyr::select(subject, group, subgroup, moca_total, age)
  
  n_Fin <- n_Fin %>%
    inner_join(Fin_Par, by="subject")
  n_Fin$subject <- paste0("f-",n_Fin$subject)
  
  
  
  #Iowa metadata
  Iow_Par <- read.table(file = "participants_iowa.tsv", sep = "\t", header = TRUE)
  Iow_Par <- Iow_Par %>%
    separate(participant_id, c("Sub", "subject"), "-")
  Iow_Par$subject <- as.numeric(Iow_Par$subject)
  
  Iow_Par <- Iow_Par %>%
    dplyr::select(subject, group, age, gender)
  #  dplyr::select(subject, group, subgroup, moca_total, age)
  
  n_Iow <- n_Iow %>%
    inner_join(Iow_Par, by="subject")
  n_Iow$subject <- paste0("i-",n_Iow$subject)
  
  
  
  
  ### CHANNELS IN EACH DATASET 32 COMMON IN ALL EXCEPT BY IOWA WITH 29 CHANNS
  
  cal <-n_Cal %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Pz", "Fp2", 
                                                              "O1", "O2", "FC2", "Fz", 
                                                              "F4", "P3", "Oz", "F7", 
                                                              "PO3", "CP6", "FC1", "P8", 
                                                              "PO4", "T7", "Cz", "F3",  
                                                              "CP2","FC6", "C3", "FC5", 
                                                              "CP5", "AF4", "F8", "P4",  
                                                              "CP1", "C4", "AF3", "T8",  
                                                              "P7")))
  
  med <-n_Med %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Pz", "Fp2", 
                                                              "O1", "O2", "FC2", "Fz", 
                                                              "F4", "P3", "Oz", "F7", 
                                                              "PO3", "CP6", "FC1", "P8", 
                                                              "PO4", "T7", "Cz", "F3",  
                                                              "CP2","FC6", "C3", "FC5", 
                                                              "CP5", "AF4", "F8", "P4",  
                                                              "CP1", "C4", "AF3", "T8",  
                                                              "P7")))
  fin <-n_Fin %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Pz", "Fp2", 
                                                              "O1", "O2", "FC2", "Fz", 
                                                              "F4", "P3", "Oz", "F7", 
                                                              "PO3", "CP6", "FC1", "P8", 
                                                              "PO4", "T7", "Cz", "F3",  
                                                              "CP2","FC6", "C3", "FC5", 
                                                              "CP5", "AF4", "F8", "P4",  
                                                              "CP1", "C4", "AF3", "T8",  
                                                              "P7")))
  
  
  
  iow <-n_Iow %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2", "O1", "O2", "FC2", "Fz", "F4", "P3", "Oz", 
                                                              "F7", "CP6", "FC1", "P8", "T7", "Cz", "F3",  "CP2", "FC6", 
                                                              "C3", "FC5", "CP5", "AF4", "F8", "P4",  "CP1", "C4", "AF3", 
                                                              "T8",  "P7")))
  
  
  ##################################################### IOWA DATAFRAME ####
  
  iow$d_log <- log10(iow$delta)
  iow$t_log <- log10(iow$theta)
  iow$st_log <- log10(iow$slow_theta)
  iow$a_log <- log10(iow$alpha)
  iow$b_log <- log10(iow$beta)
  iow$prealpha_log <- log10(iow$pre_alpha)
  iow$atr_log <- log10(iow$alpha_theta)
  
  
  mean_prealpha <- iow %>% group_by(subject, channel) %>% summarise (prealpha_mean = mean(prealpha_log))
  mean_d_log <- iow %>% group_by(subject, channel) %>% summarise (d_mean = mean(d_log))
  mean_t_log <- iow %>% group_by(subject, channel) %>% summarise (t_mean = mean(t_log))
  mean_st_log <- iow %>% group_by(subject, channel) %>% summarise (st_mean = mean(st_log))
  
  mean_a_log <- iow %>% group_by(subject, channel) %>% summarise (a_mean = mean(a_log))
  mean_b_log <- iow %>% group_by(subject, channel) %>% summarise (b_mean = mean(b_log))
  mean_atr_log <- iow %>% group_by(subject, channel) %>% summarise (atr_mean = mean(atr_log))
  
  
  median_prealpha <- iow %>% group_by(subject, channel) %>% summarise (prealpha_median = median(prealpha_log))
  median_d_log <- iow %>% group_by(subject, channel) %>% summarise (d_median = median(d_log))
  median_t_log <- iow %>% group_by(subject, channel) %>% summarise (t_median = median(t_log))
  median_st_log <- iow %>% group_by(subject, channel) %>% summarise (st_median = median(st_log))
  median_a_log <- iow %>% group_by(subject, channel) %>% summarise (a_median = median(a_log))
  median_b_log <- iow %>% group_by(subject, channel) %>% summarise (b_median = median(b_log))
  median_atr_log <- iow %>% group_by(subject, channel) %>% summarise (atr_median = median(atr_log))
  
  
  ### MEAN SPECTRAL ####
  
  #MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
  mean_spectral <- iow  %>%  filter(epoch == 0)
  
  mean_spectral <- mean_spectral %>%
    dplyr::select(subject, channel, group, gender)
  #  select(subject, group, subgroup, moca_total, age)
  
  ##################################################### BOXPLOTS IOWA ########################################
  ### MEAN SPECTRAL ####
  mean_spectral <- mean_spectral %>%
    inner_join(mean_prealpha, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(mean_d_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(mean_t_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(mean_st_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(mean_a_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(mean_b_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(mean_atr_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_atr_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_prealpha, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_a_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_t_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_st_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_d_log, by=c("subject","channel"))
  
  mean_spectral <- mean_spectral %>%
    inner_join(median_b_log, by=c("subject","channel"))

  #### TOPOPLOTS ####
  group_prealpha <- mean_spectral %>% group_by(channel, group) %>% summarise (prealpha_group = mean(prealpha_mean))
  group_slowtheta <- mean_spectral %>% group_by(channel, group) %>% summarise (s_theta_group = mean(st_mean))
  group_delta <- mean_spectral %>% group_by(channel, group) %>% summarise (delta_group = mean(d_mean))
  group_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (theta_group = mean(t_mean))
  group_alpha <- mean_spectral %>% group_by(channel, group) %>% summarise (alpha_group = mean(a_mean))
  group_beta <- mean_spectral %>% group_by(channel, group) %>% summarise (beta_group = mean(b_mean))
  group_alpha_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (mean_atr_group = mean(atr_mean))
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(group_prealpha, file ="iowa_group_prealpha.csv", row.names=FALSE)
  write.csv(group_slowtheta, file ="iowa_group_slow_theta.csv", row.names=FALSE)
  write.csv(group_delta, file ="iowa_group_delta.csv", row.names=FALSE)
  write.csv(group_theta, file ="iowa_group_theta.csv", row.names=FALSE)
  write.csv(group_alpha, file ="iowa_group_alpha.csv", row.names=FALSE)
  write.csv(group_beta, file ="iowa_group_beta.csv", row.names=FALSE)
  write.csv(group_alpha_theta, file ="iowa_group_mean_alpha_theta.csv", row.names=FALSE)
  
  

  #### MEAN DELTA - NOT SIGNIFICANT ALL CHANNELS  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(d_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "d_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Delta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_delta.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(d_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_delta.csv", row.names=FALSE)
  
  
  #### MEAN THETA - SIGNIFICANT ALL CHANNELS  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(t_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "t_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_theta.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(t_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_theta.csv", row.names=FALSE)
  
    
  #### MEAN SLOW THETA - NOT SIGNIFICANT ALL CHANNELS  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(st_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "st_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Slow-theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_slow_theta.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(st_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_slow_theta.csv", row.names=FALSE)
  
  
  #### MEAN PREALPHA - SIGNIFICANT ALL CHANNELS  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(prealpha_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "prealpha_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Pre-alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_prealpha.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(prealpha_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_prealpha.csv", row.names=FALSE)
  
  
  #### MEAN ALPHA - NOT SIGNIFICANT ALL CHANNELS  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(a_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "a_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_alpha.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(a_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_alpha.csv", row.names=FALSE)
  
  
  #### MEAN BETA - SIGNIFICANT MOST CHANNELS  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(b_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "b_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Beta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_beta.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(b_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_beta.csv", row.names=FALSE)
  
  
  #### MEAN ALPHA/THETA - NOT SIGNIFICANT ALL CHANNELS, SOME SIGNIFICANT UNCORRECTED P  ####
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(atr_mean ~ group, ref.group = "hc") %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  
  stat.test <- stat.test %>%
    mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))
  
  stat.test <- stat.test %>% add_xy_position(x = "group")
  
  
  # Create a box plot
  bxp <- ggboxplot(
    mean_spectral, x = "group", y = "atr_mean" , color = "group",  add = "jitter",
    facet.by = "channel",  xlab = F,
    ylab = "Mean Log 10 (Alpha/Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))
  
  
  
  # Make facet and add p-values
  # Gráficos
  setwd("D:/escritorio/tSNE_entropia_R/Publication")
  png(filename = "./Iowa_atr.png", width = 16, height = 10, units = "in", res = 300)
  bxp + 
    stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                       hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_x_discrete(labels = c("non-PD", "PD"))+
    scale_color_discrete(name="Group",
                         labels=c("non-PD","PD"), direction = -1) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"))
  dev.off()
  
  stat.test <- mean_spectral %>%
    group_by(channel) %>%
    t_test(atr_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance() %>%
    mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))
  
  
  stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
  stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)
  
  
  setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")
  
  write.csv(stat.test, file ="iowa_topoplot_atr.csv", row.names=FALSE)
  
  
  #### BY GENDER - MEAN PREALPHA ####

stat.test <- mean_spectral %>%
  group_by(channel, group) %>%
  t_test(prealpha_mean ~ gender , ref.group = "f") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "gender",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) PSD", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.05, hide.ns = F,
                     label = "{p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Gender",
                     labels=c("Female","Male"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))





##################################################### CALIFORNIA DATAFRAME ####

cal$d_log <- log10(cal$delta)
cal$t_log <- log10(cal$theta)
cal$st_log <- log10(cal$slow_theta)
cal$a_log <- log10(cal$alpha)
cal$b_log <- log10(cal$beta)
cal$prealpha_log <- log10(cal$pre_alpha)
cal$atr_log <- log10(cal$alpha_theta)


mean_prealpha <- cal %>% group_by(subject, channel) %>% summarise (prealpha_mean = mean(prealpha_log))
mean_d_log <- cal %>% group_by(subject, channel) %>% summarise (d_mean = mean(d_log))
mean_t_log <- cal %>% group_by(subject, channel) %>% summarise (t_mean = mean(t_log))
mean_st_log <- cal %>% group_by(subject, channel) %>% summarise (st_mean = mean(st_log))

mean_a_log <- cal %>% group_by(subject, channel) %>% summarise (a_mean = mean(a_log))
mean_b_log <- cal %>% group_by(subject, channel) %>% summarise (b_mean = mean(b_log))
mean_atr_log <- cal %>% group_by(subject, channel) %>% summarise (atr_mean = mean(atr_log))


median_prealpha <- cal %>% group_by(subject, channel) %>% summarise (prealpha_median = median(prealpha_log))
median_d_log <- cal %>% group_by(subject, channel) %>% summarise (d_median = median(d_log))
median_t_log <- cal %>% group_by(subject, channel) %>% summarise (t_median = median(t_log))
median_st_log <- cal %>% group_by(subject, channel) %>% summarise (st_median = median(st_log))
median_a_log <- cal %>% group_by(subject, channel) %>% summarise (a_median = median(a_log))
median_b_log <- cal %>% group_by(subject, channel) %>% summarise (b_median = median(b_log))
median_atr_log <- cal %>% group_by(subject, channel) %>% summarise (atr_median = median(atr_log))


### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_spectral <- cal  %>%  filter(epoch == 0)

mean_spectral <- mean_spectral %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS CALIFORNIA ########################################
### MEAN SPECTRAL ####
mean_spectral <- mean_spectral %>%
  inner_join(mean_prealpha, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_d_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_t_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_st_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_a_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_b_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_atr_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_atr_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_prealpha, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_a_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_t_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_st_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_d_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_b_log, by=c("subject","channel"))

#### TOPOPLOTS ####
group_prealpha <- mean_spectral %>% group_by(channel, group) %>% summarise (prealpha_group = mean(prealpha_mean))
group_slowtheta <- mean_spectral %>% group_by(channel, group) %>% summarise (s_theta_group = mean(st_mean))
group_delta <- mean_spectral %>% group_by(channel, group) %>% summarise (delta_group = mean(d_mean))
group_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (theta_group = mean(t_mean))
group_alpha <- mean_spectral %>% group_by(channel, group) %>% summarise (alpha_group = mean(a_mean))
group_beta <- mean_spectral %>% group_by(channel, group) %>% summarise (beta_group = mean(b_mean))
group_alpha_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (mean_atr_group = mean(atr_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_prealpha, file ="california_group_prealpha.csv", row.names=FALSE)
write.csv(group_slowtheta, file ="california_group_slow_theta.csv", row.names=FALSE)
write.csv(group_delta, file ="california_group_delta.csv", row.names=FALSE)
write.csv(group_theta, file ="california_group_theta.csv", row.names=FALSE)
write.csv(group_alpha, file ="california_group_alpha.csv", row.names=FALSE)
write.csv(group_beta, file ="california_group_beta.csv", row.names=FALSE)
write.csv(group_alpha_theta, file ="california_group_mean_alpha_theta.csv", row.names=FALSE)

#### MEAN DELTA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(d_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "d_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Delta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_delta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(d_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_delta.csv", row.names=FALSE)


#### MEAN THETA - SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(t_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "t_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_theta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(t_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_theta.csv", row.names=FALSE)


#### MEAN SLOW THETA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(st_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "st_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Slow-theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_slow_theta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(st_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_slow_theta.csv", row.names=FALSE)


#### MEAN PREALPHA - SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(prealpha_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_prealpha.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(prealpha_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_prealpha.csv", row.names=FALSE)


#### MEAN ALPHA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(a_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "a_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_alpha.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(a_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_alpha.csv", row.names=FALSE)


#### MEAN BETA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(b_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "b_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Beta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_beta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(b_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_beta.csv", row.names=FALSE)


#### MEAN ALPHA/THETA - NOT SIGNIFICANT ALL CHANNELS, NONE SIGNIFICANT UNCORRECTED P  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(atr_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "atr_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Alpha/Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_atr.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()



stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(atr_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_atr.csv", row.names=FALSE)


#### BY GENDER - MEAN PREALPHA ####

stat.test <- mean_spectral %>%
  group_by(channel, group) %>%
  t_test(prealpha_mean ~ gender , ref.group = "f") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "gender",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) PSD", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.05, hide.ns = F,
                     label = "{p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Gender",
                       labels=c("Female","Male"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))






##################################################### FINLAND DATAFRAME ####

fin$d_log <- log10(fin$delta)
fin$t_log <- log10(fin$theta)
fin$st_log <- log10(fin$slow_theta)
fin$a_log <- log10(fin$alpha)
fin$b_log <- log10(fin$beta)
fin$prealpha_log <- log10(fin$pre_alpha)
fin$atr_log <- log10(fin$alpha_theta)


mean_prealpha <- fin %>% group_by(subject, channel) %>% summarise (prealpha_mean = mean(prealpha_log))
mean_d_log <- fin %>% group_by(subject, channel) %>% summarise (d_mean = mean(d_log))
mean_t_log <- fin %>% group_by(subject, channel) %>% summarise (t_mean = mean(t_log))
mean_st_log <- fin %>% group_by(subject, channel) %>% summarise (st_mean = mean(st_log))

mean_a_log <- fin %>% group_by(subject, channel) %>% summarise (a_mean = mean(a_log))
mean_b_log <- fin %>% group_by(subject, channel) %>% summarise (b_mean = mean(b_log))
mean_atr_log <- fin %>% group_by(subject, channel) %>% summarise (atr_mean = mean(atr_log))


median_prealpha <- fin %>% group_by(subject, channel) %>% summarise (prealpha_median = median(prealpha_log))
median_d_log <- fin %>% group_by(subject, channel) %>% summarise (d_median = median(d_log))
median_t_log <- fin %>% group_by(subject, channel) %>% summarise (t_median = median(t_log))
median_st_log <- fin %>% group_by(subject, channel) %>% summarise (st_median = median(st_log))
median_a_log <- fin %>% group_by(subject, channel) %>% summarise (a_median = median(a_log))
median_b_log <- fin %>% group_by(subject, channel) %>% summarise (b_median = median(b_log))
median_atr_log <- fin %>% group_by(subject, channel) %>% summarise (atr_median = median(atr_log))


### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_spectral <- fin  %>%  filter(epoch == 0)

mean_spectral <- mean_spectral %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS FINLAND ########################################
### MEAN SPECTRAL ####
mean_spectral <- mean_spectral %>%
  inner_join(mean_prealpha, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_d_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_t_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_st_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_a_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_b_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_atr_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_atr_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_prealpha, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_a_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_t_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_st_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_d_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_b_log, by=c("subject","channel"))

#### TOPOPLOTS ####
group_prealpha <- mean_spectral %>% group_by(channel, group) %>% summarise (prealpha_group = mean(prealpha_mean))
group_slowtheta <- mean_spectral %>% group_by(channel, group) %>% summarise (s_theta_group = mean(st_mean))
group_delta <- mean_spectral %>% group_by(channel, group) %>% summarise (delta_group = mean(d_mean))
group_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (theta_group = mean(t_mean))
group_alpha <- mean_spectral %>% group_by(channel, group) %>% summarise (alpha_group = mean(a_mean))
group_beta <- mean_spectral %>% group_by(channel, group) %>% summarise (beta_group = mean(b_mean))
group_alpha_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (mean_atr_group = mean(atr_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_prealpha, file ="finland_group_prealpha.csv", row.names=FALSE)
write.csv(group_slowtheta, file ="finland_group_slow_theta.csv", row.names=FALSE)
write.csv(group_delta, file ="finland_group_delta.csv", row.names=FALSE)
write.csv(group_theta, file ="finland_group_theta.csv", row.names=FALSE)
write.csv(group_alpha, file ="finland_group_alpha.csv", row.names=FALSE)
write.csv(group_beta, file ="finland_group_beta.csv", row.names=FALSE)
write.csv(group_alpha_theta, file ="finland_group_mean_alpha_theta.csv", row.names=FALSE)

#### MEAN DELTA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(d_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "d_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Delta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_delta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(d_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_delta.csv", row.names=FALSE)


#### MEAN THETA - NOT SO SIGNIFICANT ONLY SOME CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(t_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "t_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_theta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(t_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_theta.csv", row.names=FALSE)


#### MEAN SLOW THETA - NOT SIGNIFICANT MOST CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(st_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "st_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Slow-theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_slow_theta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(st_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_slow_theta.csv", row.names=FALSE)


#### MEAN PREALPHA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(prealpha_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_prealpha.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(prealpha_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_prealpha.csv", row.names=FALSE)


#### MEAN ALPHA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(a_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "a_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_alpha.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(a_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_alpha.csv", row.names=FALSE)


#### MEAN BETA - NOT SIGNIFICANT MOST CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(b_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "b_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Beta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_beta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(b_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_beta.csv", row.names=FALSE)


#### MEAN ALPHA/THETA - MOST SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(atr_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "atr_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Alpha/Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_atr.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(atr_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_atr.csv", row.names=FALSE)


#### BY GENDER - MEAN PREALPHA ####

stat.test <- mean_spectral %>%
  group_by(channel, group) %>%
  t_test(prealpha_mean ~ gender , ref.group = "f") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "gender",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) PSD", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.05, hide.ns = F,
                     label = "{p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Gender",
                       labels=c("Female","Male"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))







##################################################### MEDELLIN DATAFRAME ####

med$d_log <- log10(med$delta)
med$t_log <- log10(med$theta)
med$st_log <- log10(med$slow_theta)
med$a_log <- log10(med$alpha)
med$b_log <- log10(med$beta)
med$prealpha_log <- log10(med$pre_alpha)
med$atr_log <- log10(med$alpha_theta)


mean_prealpha <- med %>% group_by(subject, channel) %>% summarise (prealpha_mean = mean(prealpha_log))
mean_d_log <- med %>% group_by(subject, channel) %>% summarise (d_mean = mean(d_log))
mean_t_log <- med %>% group_by(subject, channel) %>% summarise (t_mean = mean(t_log))
mean_st_log <- med %>% group_by(subject, channel) %>% summarise (st_mean = mean(st_log))

mean_a_log <- med %>% group_by(subject, channel) %>% summarise (a_mean = mean(a_log))
mean_b_log <- med %>% group_by(subject, channel) %>% summarise (b_mean = mean(b_log))
mean_atr_log <- med %>% group_by(subject, channel) %>% summarise (atr_mean = mean(atr_log))


median_prealpha <- med %>% group_by(subject, channel) %>% summarise (prealpha_median = median(prealpha_log))
median_d_log <- med %>% group_by(subject, channel) %>% summarise (d_median = median(d_log))
median_t_log <- med %>% group_by(subject, channel) %>% summarise (t_median = median(t_log))
median_st_log <- med %>% group_by(subject, channel) %>% summarise (st_median = median(st_log))
median_a_log <- med %>% group_by(subject, channel) %>% summarise (a_median = median(a_log))
median_b_log <- med %>% group_by(subject, channel) %>% summarise (b_median = median(b_log))
median_atr_log <- med %>% group_by(subject, channel) %>% summarise (atr_median = median(atr_log))


### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_spectral <- med  %>%  filter(epoch == 0)

mean_spectral <- mean_spectral %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS MEDELLIN ########################################
### MEAN SPECTRAL ####
mean_spectral <- mean_spectral %>%
  inner_join(mean_prealpha, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_d_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_t_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_st_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_a_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_b_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(mean_atr_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_atr_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_prealpha, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_a_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_t_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_st_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_d_log, by=c("subject","channel"))

mean_spectral <- mean_spectral %>%
  inner_join(median_b_log, by=c("subject","channel"))

#### TOPOPLOTS ####
group_prealpha <- mean_spectral %>% group_by(channel, group) %>% summarise (prealpha_group = mean(prealpha_mean))
group_slowtheta <- mean_spectral %>% group_by(channel, group) %>% summarise (s_theta_group = mean(st_mean))
group_delta <- mean_spectral %>% group_by(channel, group) %>% summarise (delta_group = mean(d_mean))
group_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (theta_group = mean(t_mean))
group_alpha <- mean_spectral %>% group_by(channel, group) %>% summarise (alpha_group = mean(a_mean))
group_beta <- mean_spectral %>% group_by(channel, group) %>% summarise (beta_group = mean(b_mean))
group_alpha_theta <- mean_spectral %>% group_by(channel, group) %>% summarise (mean_atr_group = mean(atr_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_prealpha, file ="medellin_group_prealpha.csv", row.names=FALSE)
write.csv(group_slowtheta, file ="medellin_group_slow_theta.csv", row.names=FALSE)
write.csv(group_delta, file ="medellin_group_delta.csv", row.names=FALSE)
write.csv(group_theta, file ="medellin_group_theta.csv", row.names=FALSE)
write.csv(group_alpha, file ="medellin_group_alpha.csv", row.names=FALSE)
write.csv(group_beta, file ="medellin_group_beta.csv", row.names=FALSE)
write.csv(group_alpha_theta, file ="medellin_group_mean_alpha_theta.csv", row.names=FALSE)


#### MEAN DELTA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(d_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "d_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Delta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_delta.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(d_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_delta.csv", row.names=FALSE)


#### MEAN THETA - SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(t_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "t_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_theta.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(t_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_theta.csv", row.names=FALSE)


#### MEAN SLOW THETA -  SIGNIFICANT SOME CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(st_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "st_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Slow-theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_slow_theta.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(st_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_slow_theta.csv", row.names=FALSE)


#### MEAN PREALPHA - SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(prealpha_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_prealpha.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(prealpha_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_prealpha.csv", row.names=FALSE)


#### MEAN ALPHA - NOT SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(a_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "a_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Alpha) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_alpha.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(a_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_alpha.csv", row.names=FALSE)


#### MEAN BETA - NOT SIGNIFICANT MOST CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(b_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "b_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Beta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_beta.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(b_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_beta.csv", row.names=FALSE)


#### MEAN ALPHA/THETA -  SIGNIFICANT ALL CHANNELS  ####

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(atr_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "atr_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Alpha/Theta) relative PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))



# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_atr.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()

stat.test <- mean_spectral %>%
  group_by(channel) %>%
  t_test(atr_mean ~ group, ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))


stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_atr.csv", row.names=FALSE)


#### BY GENDER - MEAN PREALPHA ####

stat.test <- mean_spectral %>%
  group_by(channel, group) %>%
  t_test(prealpha_mean ~ gender , ref.group = "f") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_spectral, x = "group", y = "prealpha_mean" , color = "gender",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Log 10 (Pre-alpha) PSD", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.05, hide.ns = F,
                     label = "{p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Gender",
                       labels=c("Female","Male"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))








###############################################################################################   DOMINANT FREQ ########
#### READ CSVs WITH SPECTRAL FEATURES ####
setwd(dir = "D:/escritorio/tSNE_entropia_R")

n_Cal <- read.table(file = "df_california_5s_no_overlapping.csv", sep = ",", header = TRUE)
n_Med <- read.table(file = "df_medellin_5s_no_overlapping.csv", sep = ",", header = TRUE)
n_Fin <- read.table(file = "df_finland_5s_no_overlapping.csv", sep = ",", header = TRUE)
n_Iow <- read.table(file = "df_iowa_5s_no_overlapping.csv", sep = ",", header = TRUE)




#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
Med_Par <- read.table(file = "participants_med.tsv", sep = "\t", header = TRUE)
Med_Par <- Med_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Med_Par$subject <- as.numeric(Med_Par$subject)

Med_Par <- Med_Par %>%
  dplyr::select(subject, group, age, gender)
#  dplyr::select(subject, group, subgroup, moca_total, age)

n_Med <- n_Med %>%
  inner_join(Med_Par, by="subject")
n_Med$subject <- paste0("m-",n_Med$subject)



#California Metadata
Cal_par <- read.table(file = "participants_cal.tsv", sep = "\t", header = TRUE)
Cal_par <- Cal_par %>%
  separate(participant_id, c("Sub", "subject"), "-")

Cal_par <- Cal_par %>%
  dplyr::select(subject, age, gender, group)
#  dplyr::select(subject, group, subgroup, moca_total, age)

n_Cal <- n_Cal %>%
  inner_join(Cal_par, by="subject")

n_Cal$subject <- paste0("c-",n_Cal$subject)




#Finland metadata
Fin_Par <- read.table(file = "participants_fin.tsv", sep = "\t", header = TRUE)
Fin_Par <- Fin_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Fin_Par$subject <- as.numeric(Fin_Par$subject)

Fin_Par <- Fin_Par %>%
  dplyr::select(subject, group, age, gender)
#  dplyr::select(subject, group, subgroup, moca_total, age)

n_Fin <- n_Fin %>%
  inner_join(Fin_Par, by="subject")
n_Fin$subject <- paste0("f-",n_Fin$subject)



#Iowa metadata
Iow_Par <- read.table(file = "participants_iowa.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  dplyr::select(subject, group, subgroup, moca_total, age)

n_Iow <- n_Iow %>%
  inner_join(Iow_Par, by="subject")
n_Iow$subject <- paste0("i-",n_Iow$subject)




### CHANNELS IN EACH DATASET 32 COMMON IN ALL EXCEPT BY IOWA WITH 29 CHANNS

cal <-n_Cal %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Pz", "Fp2", 
                                                            "O1", "O2", "FC2", "Fz", 
                                                            "F4", "P3", "Oz", "F7", 
                                                            "PO3", "CP6", "FC1", "P8", 
                                                            "PO4", "T7", "Cz", "F3",  
                                                            "CP2","FC6", "C3", "FC5", 
                                                            "CP5", "AF4", "F8", "P4",  
                                                            "CP1", "C4", "AF3", "T8",  
                                                            "P7")))

med <-n_Med %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Pz", "Fp2", 
                                                            "O1", "O2", "FC2", "Fz", 
                                                            "F4", "P3", "Oz", "F7", 
                                                            "PO3", "CP6", "FC1", "P8", 
                                                            "PO4", "T7", "Cz", "F3",  
                                                            "CP2","FC6", "C3", "FC5", 
                                                            "CP5", "AF4", "F8", "P4",  
                                                            "CP1", "C4", "AF3", "T8",  
                                                            "P7")))
fin <-n_Fin %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Pz", "Fp2", 
                                                            "O1", "O2", "FC2", "Fz", 
                                                            "F4", "P3", "Oz", "F7", 
                                                            "PO3", "CP6", "FC1", "P8", 
                                                            "PO4", "T7", "Cz", "F3",  
                                                            "CP2","FC6", "C3", "FC5", 
                                                            "CP5", "AF4", "F8", "P4",  
                                                            "CP1", "C4", "AF3", "T8",  
                                                            "P7")))



iow <-n_Iow %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2", "O1", "O2", "FC2", "Fz", "F4", "P3", "Oz", 
                                                            "F7", "CP6", "FC1", "P8", "T7", "Cz", "F3",  "CP2", "FC6", 
                                                            "C3", "FC5", "CP5", "AF4", "F8", "P4",  "CP1", "C4", "AF3", 
                                                            "T8",  "P7")))


##################################################### IOWA DATAFRAME ####


mean_df <- iow %>% group_by(subject, channel) %>% summarise (df_mean = mean(df_ep))
mean_dfv <- iow %>% group_by(subject, channel) %>% summarise (dfv_mean = mean(dfv_ch))
median_df <- iow %>% group_by(subject, channel) %>% summarise (df_median = median(df_ep))
median_dfv <- iow %>% group_by(subject, channel) %>% summarise (dfv_median = median(dfv_ch))



### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_dominant <- iow  %>%  filter(epoch == 0)

mean_dominant <- mean_dominant %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS IOWA ########################################
### MEAN SPECTRAL ####
mean_dominant <- mean_dominant %>%
  inner_join(mean_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(mean_dfv, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_dfv, by=c("subject","channel"))


#### TOPOPLOTS ####
group_df <- mean_dominant %>% group_by(channel, group) %>% summarise (df_group = mean(df_mean))
group_dfv <- mean_dominant %>% group_by(channel, group) %>% summarise (dfv_group = mean(dfv_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_df, file ="iowa_group_df.csv", row.names=FALSE)
write.csv(group_dfv, file ="iowa_group_dfv.csv", row.names=FALSE)


#### MEAN DF CHANNELS SOME RIGHT PARIETAL SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Iowa_df.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="iowa_topoplot_df.csv", row.names=FALSE)

#### MEAN DFV SOME SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Iowa_dfv.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="iowa_topoplot_dfv.csv", row.names=FALSE)






#### MEDIAN DF CHANNELS NONE SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Iowa", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Iowa_df_median.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()



##################################################### CALIFORNIA DATAFRAME ####


mean_df <- cal %>% group_by(subject, channel) %>% summarise (df_mean = mean(df_ep))
mean_dfv <- cal %>% group_by(subject, channel) %>% summarise (dfv_mean = mean(dfv_ch))
median_df <- cal %>% group_by(subject, channel) %>% summarise (df_median = median(df_ep))
median_dfv <- cal %>% group_by(subject, channel) %>% summarise (dfv_median = median(dfv_ch))



### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_dominant <- cal  %>%  filter(epoch == 0)

mean_dominant <- mean_dominant %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS CALIFORNIA ########################################
### MEAN SPECTRAL ####
mean_dominant <- mean_dominant %>%
  inner_join(mean_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(mean_dfv, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_dfv, by=c("subject","channel"))


#### TOPOPLOTS ####
group_df <- mean_dominant %>% group_by(channel, group) %>% summarise (df_group = mean(df_mean))
group_dfv <- mean_dominant %>% group_by(channel, group) %>% summarise (dfv_group = mean(dfv_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_df, file ="california_group_df.csv", row.names=FALSE)
write.csv(group_dfv, file ="california_group_dfv.csv", row.names=FALSE)


#### MEAN DF CHANNELS NONE SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_df.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_df.csv", row.names=FALSE)

#### MEAN DFV NONE SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_dfv.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="california_topoplot_dfv.csv", row.names=FALSE)





#### MEDIAN DF CHANNELS SOME SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_df_median.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} p = {p} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


#### MEDIAN DFV CHANNELS SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "California", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./California_dfv_median.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()











##################################################### FINLAND DATAFRAME ####
setwd("D:/escritorio/tSNE_entropia_R")

mean_df <- fin %>% group_by(subject, channel) %>% summarise (df_mean = mean(df_ep))
mean_dfv <- fin %>% group_by(subject, channel) %>% summarise (dfv_mean = mean(dfv_ch))
median_df <- fin %>% group_by(subject, channel) %>% summarise (df_median = median(df_ep))
median_dfv <- fin %>% group_by(subject, channel) %>% summarise (dfv_median = median(dfv_ch))




### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_dominant <- fin  %>%  filter(epoch == 0)

mean_dominant <- mean_dominant %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS FINLAND ########################################
### MEAN SPECTRAL ####
mean_dominant <- mean_dominant %>%
  inner_join(mean_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(mean_dfv, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_dfv, by=c("subject","channel"))

#### TOPOPLOTS ####
group_df <- mean_dominant %>% group_by(channel, group) %>% summarise (df_group = mean(df_mean))
group_dfv <- mean_dominant %>% group_by(channel, group) %>% summarise (dfv_group = mean(dfv_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_df, file ="finland_group_df.csv", row.names=FALSE)
write.csv(group_dfv, file ="finland_group_dfv.csv", row.names=FALSE)

#### MEAN DF CHANNELS NONE SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_df.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_df.csv", row.names=FALSE)

#### MEAN DFV NONE SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_dfv.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="finland_topoplot_dfv.csv", row.names=FALSE)








#### MEDIAN DF CHANNELS SOME SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_df_median.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} p = {p} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


#### MEDIAN DFV CHANNELS SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Finland", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Finland_dfv_median.png", width = 16, height = 10, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()




##################################################### MEDELLIN DATAFRAME ####
setwd("D:/escritorio/tSNE_entropia_R ")


mean_df <- med %>% group_by(subject, channel) %>% summarise (df_mean = mean(df_ep))
mean_dfv <- med %>% group_by(subject, channel) %>% summarise (dfv_mean = mean(dfv_ch))
median_df <- med %>% group_by(subject, channel) %>% summarise (df_median = median(df_ep))
median_dfv <- med %>% group_by(subject, channel) %>% summarise (dfv_median = median(dfv_ch))



### MEAN SPECTRAL ####

#MERGE PARTICIPANTS METADATA: GROUP AND AGE VARIABLES
mean_dominant <- med  %>%  filter(epoch == 0)

mean_dominant <- mean_dominant %>%
  dplyr::select(subject, channel, group, gender)
#  select(subject, group, subgroup, moca_total, age)

##################################################### BOXPLOTS MEDELLIN ########################################
### MEAN SPECTRAL ####
mean_dominant <- mean_dominant %>%
  inner_join(mean_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(mean_dfv, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_df, by=c("subject","channel"))

mean_dominant <- mean_dominant %>%
  inner_join(median_dfv, by=c("subject","channel"))

#### TOPOPLOTS ####
group_df <- mean_dominant %>% group_by(channel, group) %>% summarise (df_group = mean(df_mean))
group_dfv <- mean_dominant %>% group_by(channel, group) %>% summarise (dfv_group = mean(dfv_mean))


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(group_df, file ="medellin_group_df.csv", row.names=FALSE)
write.csv(group_dfv, file ="medellin_group_dfv.csv", row.names=FALSE)

#### MEAN DF CHANNELS SOME SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_df.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_df.csv", row.names=FALSE)

#### MEAN DFV NONE SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_mean" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Mean Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_dfv.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_mean ~ group,  ref.group = "hc", comparisons = c("pd", "hc"), detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj,  digits = 2, nsmall = 1)))

stat.test$sign <- ifelse(stat.test$p.adj.signif < 0.05, 1, 0)
stat.test$sign_unc <- ifelse(stat.test$p < 0.05, 1, 0)


setwd(dir = "D:/escritorio/tSNE_entropia_R/TOPOPLOT/data/")

write.csv(stat.test, file ="medellin_topoplot_dfv.csv", row.names=FALSE)











#### MEDIAN DF CHANNELS SOME SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(df_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "df_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_df_median.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()


#### MEDIAN DFV CHANNELS SIGNIFICANT ####

stat.test <- mean_dominant %>%
  group_by(channel) %>%
  t_test(dfv_median ~ group, ref.group = "hc") %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = F, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  mean_dominant, x = "group", y = "dfv_median" , color = "group",  add = "jitter",
  facet.by = "channel",  xlab = F,
  ylab = "Median Dominant Frequency Variability (Hz)", title = "PD vs. Non-PD channel differences", subtitle = "Medellin", ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
# Gráficos
setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "./Medellin_dfv_median.png", width = 16, height = 11, units = "in", res = 300)
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.x = -2, 
                     hide.ns = F, label = "{labels.T} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Group",
                       labels=c("non-PD","PD"), direction = -1) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))
dev.off()












################# LEDD, UPDRS, AGE BY CENTER ######################
setwd("D:/escritorio/tSNE_entropia_R/Publication")
base <- read.table(file = "clinical_dataset_all.csv", sep = ";", header = TRUE)
stat.test <- base %>%
  group_by(group) %>%
  t_test(age ~ center) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "group")


# Create a box plot
bxp <- ggboxplot(
  base, x = "group", y = "age" , color = "center",  add = "jitter", xlab = F,
  ylab = "Age (years)", title = "Age in PD and Non-PD subjects by center",  ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.02, hide.ns = F,
                     label = "p = {p} FDR = {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(labels = c("non-PD", "PD"))+
  scale_color_discrete(name="Dataset",
                       labels=c("California","Finland", "Iowa","Medellin")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))





base <- base  %>%  filter(group == "pd")
sort <- base[order(base$center),]

base = sort

stat.test <- base %>%
  t_test(PD_duration ~ center) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "center")


# Create a box plot
bxp <- ggboxplot(
  base, x = "center", y = "PD_duration" , color = "center",  add = "jitter", xlab = F,
  ylab = "Duration of PD (months)", title = "Duration of PD by center",  ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.02, hide.ns = F,
                     label = "p = {p} FDR = {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_discrete(name="Dataset",
                       labels=c("California","Finland", "Iowa","Medellin")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))



stat.test <- base %>%
  t_test(updrs_iii ~ center) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "center")


# Create a box plot
bxp <- ggboxplot(
  base, x = "center", y = "updrs_iii" , color = "center",  add = "jitter", xlab = F,
  ylab = "UPDRS - III", title = "UPDRS - III in PD subjects by center",  ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.02, hide.ns = F,
                     label = "p = {p} {labels.FDR} {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_discrete(name="Dataset",
                       labels=c("California","Finland", "Iowa","Medellin")) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))



base <- base  %>%  filter(center != "Iowa")

stat.test <- base %>%
  t_test(LEDD ~ center) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(labels.FDR = paste0("FDR = ", format(p.adj, scientific = TRUE, digits = 2)))



stat.test <- stat.test %>%
  mutate(labels.T = paste0("T = ", format(stat.test$statistic, digits = 2)))

stat.test <- stat.test %>% add_xy_position(x = "center")


# Create a box plot
bxp <- ggboxplot(
  base, x = "center", y = "LEDD" , color = "center",  add = "jitter", xlab = F,
  ylab = "LEDD (mg/day)", title = "LEDD in PD subjects by center",  ggtheme = theme_pubr(border = T, base_family = "mono"))


# Make facet and add p-values
bxp + 
  stat_pvalue_manual(data = stat.test, bracket.nudge.y = 0.02, hide.ns = F,
                     label = "p = {p} FDR = {p.adj.signif}", family = "mono") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_discrete(name="Dataset")+
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"))




