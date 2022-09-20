
#Cargar paquetes
library(dplyr)
library(fda)
library(fda.usc)
library(ggplot2)
library(data.table) #Handle data.table objects
library(dplyr)      #The dplyr package provides a fairly agile way to handle R data files.
library(ggplot2)    # Graphs wuth ggplot2 library
library(reshape2)
library(Hmisc)
library(Metrics)
library(openxlsx)
require(tidyverse)
require(readxl)
library(ggpubr)

#### FIRST, ONE SHOULD EDIT TPERM.FD FUNCTION AND ADD THIS VALS TO GET NON-ABSOLUTE VALUES OF T VALS.
#fix(tperm.fd)
#Tvalores = (mean1 - mean2)/sqrt(var1 + var2)

#####################   IOWA ############################

####### RUN FROM HERE #################

#Directorio de trabajo
#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/spectral")


carpeta <- "iowa"
archivo <- "nonorm_spectral_iowa_5s_no_overlapping.csv"
indicador <- c("delta","theta","alpha"      
               ,"beta","pre_alpha","alpha_theta")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_iowa.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <- archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" ,"P4",  "P8" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  ECM_hc <- NULL
  ECM_pd <- NULL
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    ECM_hc[j] <- sqrt(mean(Sesgo_hc^2))
    ECM_pd[j] <- sqrt(mean(Sesgo_pd^2))
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  pdf(paste0("ADF_iowa/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  
  pdf(paste0("ADF_iowa/ECM_promedio",indicador[w], ".pdf"))
  plot(k,  ECM_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "ECM",
       ylim = c(min(ECM_hc, ECM_pd), max(ECM_hc, ECM_pd) )) 
  lines(k, ECM_pd, type = "b", pch = 19, col = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador
#######################################################################################
###############################################################################################################################
################################# FUNCTIONAL DATA ANALYSIS ################################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 12, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 12, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_iowa/Grafico_CurvasObs_Vs_CurvasSua_Iowa_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - Iowa",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - Iowa",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_iowa/Grafico_CurvasObs_Vs_CurvasSua_Iowa_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - Iowa - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - Iowa - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_iowa/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - Iowa",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - Iowa",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    max(TiempoEval)
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_iowa/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",29),1:29)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
}



nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8",  "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",29),1:29))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")



###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales_iowa.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "iowa_alpha_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_iowa_alpha_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/alpha.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nAlpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/iowa_beta_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_iowa_beta_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/beta.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nBeta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_higuchi_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/iowa_alpha_theta_channel_epoch.xlsx")
California_higuchi_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_iowa_alpha_theta_channel_epoch.xlsx")

California_higuchi_fd_p_long <-  California_higuchi_fd_p %>%
  mutate(Epoch = row.names(California_higuchi_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd_t_long <-  California_higuchi_fd_t %>%
  mutate(Epoch = row.names(California_higuchi_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd <- California_higuchi_fd_t_long %>%
  inner_join(y = California_higuchi_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/alpha_theta.png", width = 7, height = 5, units = "in", res = 300)
California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nAlpha/theta ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()





California_hjort_complexity_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/iowa_delta_channel_epoch.xlsx")
California_hjort_complexity_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_iowa_delta_channel_epoch.xlsx")

California_hjort_complexity_p_long <-  California_hjort_complexity_p %>%
  mutate(Epoch = row.names(California_hjort_complexity_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity_t_long <-  California_hjort_complexity_t %>%
  mutate(Epoch = row.names(California_hjort_complexity_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity <- California_hjort_complexity_t_long %>%
  inner_join(y = California_hjort_complexity_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/delta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nDelta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_hjort_mobility_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/iowa_theta_channel_epoch.xlsx")
California_hjort_mobility_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_iowa_theta_channel_epoch.xlsx")

California_hjort_mobility_p_long <-  California_hjort_mobility_p %>%
  mutate(Epoch = row.names(California_hjort_mobility_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility_t_long <-  California_hjort_mobility_t %>%
  mutate(Epoch = row.names(California_hjort_mobility_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility <- California_hjort_mobility_t_long %>%
  inner_join(y = California_hjort_mobility_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/theta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nTheta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()







California_katz_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/iowa_pre_alpha_channel_epoch.xlsx")
California_katz_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_iowa_pre_alpha_channel_epoch.xlsx")

California_katz_fd_p_long <-  California_katz_fd_p %>%
  mutate(Epoch = row.names(California_katz_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd_t_long <-  California_katz_fd_t %>%
  mutate(Epoch = row.names(California_katz_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd <- California_katz_fd_t_long %>%
  inner_join(y = California_katz_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/prealpha.png", width = 7, height = 5, units = "in", res = 300)
California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nPre-alpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()


California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nPre-alpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::turbo(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nPre-alpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  #  scale_fill_gradientn(name = "t value", colors = viridis::cubehelix(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(-5.73,5.73)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig1 = California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Delta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)




fig2 = California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Theta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)



fig3 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig4 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Beta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_bands_iowa   <- ggarrange(fig1, fig2,fig3, fig4, common.legend = T,
                                 legend = "right",
                                 ncol = 4, nrow = 1, 
                                 font.label = list(face = "bold", family = "mono"))

figure_bands_iowa   <- annotate_figure(figure_bands_iowa, top = text_grob("Iowa dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/bands_iowa.png", width = 20, height = 5, units = "in", res = 300)
figure_bands_iowa
dev.off()

fig5 = California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Pre-alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)

fig6 = California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha/theta relative PSD ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


#####################   DOMINANT IOWA  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "Iowa"
archivo <- "df_iowa_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_iowa.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <- archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" ,"P4",  "P8" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  setwd("D:/escritorio/tSNE_entropia_R/spectral/")
  
  pdf(paste0("ADF_iowa/Sesgo_promedio_dom",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador
#######################################################################################
###############################################################################################################################
################################# FUNCTIONAL DATA ANALYSIS ################################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 10, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 10, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_iowa/Grafico_CurvasObs_Vs_CurvasSua_Iowa_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - Iowa",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - Iowa",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_iowa/Grafico_CurvasObs_Vs_CurvasSua_Iowa_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - Iowa - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - Iowa - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_iowa/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - Iowa",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - Iowa",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    max(TiempoEval)
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_iowa/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",29),1:29)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8",  "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",29),1:29))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")










###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales_iowa.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "Iowa_dfv_ch_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_Iowa_dfv_ch_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/dfv.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nDominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/Iowa_df_ep_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_Iowa_df_ep_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/df.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Iowa\nDominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




fig8 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig7 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_spectral_iowa   <- ggarrange(fig6, fig5,fig7, fig8, common.legend = T,
                                    legend = "right",
                                    ncol = 4, nrow = 1, 
                                    font.label = list(face = "bold", family = "mono"))

figure_spectral_iowa   <- annotate_figure(figure_spectral_iowa, top = text_grob("Iowa dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_iowa")
png(filename = "./HEATMAPS/spectral.png", width = 20, height = 5, units = "in", res = 300)
figure_spectral_iowa
dev.off()









###############################################################################################  FINLAND ##########################
###############################################################################



#Directorio de trabajo
#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/spectral")


carpeta <- "finland"
archivo <- "nonorm_spectral_finland_5s_no_overlapping.csv"
indicador <- c("delta","theta","alpha"      
               ,"beta","pre_alpha","alpha_theta")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


# metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_fin.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  pdf(paste0("ADF_finland/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}
#termina la corrida para cada indicador#######################################################################################



##### DATOS FUNCIONALES PARA CADA MÉTRICA ####################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    #### PARA FINLANDIA EL NÚMERO ÓPTIMO DE BASES ESTÁ ENTRE 11 Y 13, SE ELIGIÓ 12
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 10, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 10, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_finland/Grafico_CurvasObs_Vs_CurvasSua_finland_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - Finland",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - Finland",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_finland/Grafico_CurvasObs_Vs_CurvasSua_finland_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - Finland - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - Finland - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_finland/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - Finland",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - Finland",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    max(TiempoEval)
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_finland/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",32),1:32)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
  
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8", "PO3", "PO4", "Pz", "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",32),1:32))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")



###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "finland_alpha_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_finland_alpha_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/alpha.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nAlpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/finland_beta_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_finland_beta_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/beta.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nBeta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_higuchi_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/finland_alpha_theta_channel_epoch.xlsx")
California_higuchi_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_finland_alpha_theta_channel_epoch.xlsx")

California_higuchi_fd_p_long <-  California_higuchi_fd_p %>%
  mutate(Epoch = row.names(California_higuchi_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd_t_long <-  California_higuchi_fd_t %>%
  mutate(Epoch = row.names(California_higuchi_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd <- California_higuchi_fd_t_long %>%
  inner_join(y = California_higuchi_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/alpha_theta.png", width = 7, height = 5, units = "in", res = 300)
California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nAlpha/theta ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()





California_hjort_complexity_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/finland_delta_channel_epoch.xlsx")
California_hjort_complexity_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_finland_delta_channel_epoch.xlsx")

California_hjort_complexity_p_long <-  California_hjort_complexity_p %>%
  mutate(Epoch = row.names(California_hjort_complexity_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity_t_long <-  California_hjort_complexity_t %>%
  mutate(Epoch = row.names(California_hjort_complexity_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity <- California_hjort_complexity_t_long %>%
  inner_join(y = California_hjort_complexity_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/delta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nDelta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_hjort_mobility_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/finland_theta_channel_epoch.xlsx")
California_hjort_mobility_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_finland_theta_channel_epoch.xlsx")

California_hjort_mobility_p_long <-  California_hjort_mobility_p %>%
  mutate(Epoch = row.names(California_hjort_mobility_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility_t_long <-  California_hjort_mobility_t %>%
  mutate(Epoch = row.names(California_hjort_mobility_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility <- California_hjort_mobility_t_long %>%
  inner_join(y = California_hjort_mobility_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/theta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nTheta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()







California_katz_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/finland_pre_alpha_channel_epoch.xlsx")
California_katz_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_finland_pre_alpha_channel_epoch.xlsx")

California_katz_fd_p_long <-  California_katz_fd_p %>%
  mutate(Epoch = row.names(California_katz_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd_t_long <-  California_katz_fd_t %>%
  mutate(Epoch = row.names(California_katz_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd <- California_katz_fd_t_long %>%
  inner_join(y = California_katz_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/prealpha.png", width = 7, height = 5, units = "in", res = 300)
California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nPre-alpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()






fig1 = California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Delta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)




fig2 = California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Theta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)



fig3 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig4 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Beta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_bands_finland   <- ggarrange(fig1, fig2,fig3, fig4, common.legend = T,
                                    legend = "right",
                                    ncol = 4, nrow = 1, 
                                    font.label = list(face = "bold", family = "mono"))

figure_bands_finland   <- annotate_figure(figure_bands_finland, top = text_grob("Finland dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/bands_finland.png", width = 20, height = 5, units = "in", res = 300)
figure_bands_finland
dev.off()

fig5 = California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Pre-alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)

fig6 = California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha/theta relative PSD ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)



#####################   DOMINANT FINLAND  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "finland"
archivo <- "df_finland_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_fin.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  setwd("D:/escritorio/tSNE_entropia_R/spectral/")
  
  pdf(paste0("ADF_finland/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador


#######################################################################################
###############################################################################################################################
################################# FUNCTIONAL DATA ANALYSIS ################################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 10, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 10, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_finland/Grafico_CurvasObs_Vs_CurvasSua_finland_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - finland",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - finland",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_finland/Grafico_CurvasObs_Vs_CurvasSua_finland_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - finland - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - finland - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_finland/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - finland",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - finland",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    max(TiempoEval)
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_finland/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",32),1:32)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8", "PO3", "PO4", "Pz", "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",32),1:32))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")










###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "finland_dfv_ch_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_finland_dfv_ch_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/dfv.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nDominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/finland_df_ep_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_finland_df_ep_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/df.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Finland\nDominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




fig8 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig7 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_spectral_finland   <- ggarrange(fig6, fig5,fig7, fig8, common.legend = T,
                                       legend = "right",
                                       ncol = 4, nrow = 1, 
                                       font.label = list(face = "bold", family = "mono"))

figure_spectral_finland   <- annotate_figure(figure_spectral_finland, top = text_grob("Finland dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_finland")
png(filename = "./HEATMAPS/spectral.png", width = 20, height = 5, units = "in", res = 300)
figure_spectral_finland
dev.off()






















###############################################################################################  MEDELLIN  ##########################
###############################################################################



#Directorio de trabajo
#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/spectral")


carpeta <- "medellin"
archivo <- "nonorm_spectral_antioquia_5s_no_overlapping.csv"
indicador <- c("delta","theta","alpha"      
               ,"beta","pre_alpha","alpha_theta")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


# metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_med.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  pdf(paste0("ADF_medellin/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}
#termina la corrida para cada indicador#######################################################################################
##### DATOS FUNCIONALES PARA CADA MÉTRICA ####################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    #### PARA MEDELLIN EL NÚMERO ÓPTIMO DE BASES ELEGIDO ES 25
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_medellin/Grafico_CurvasObs_Vs_CurvasSua_medellin_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - Medellin",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - Medellin",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_medellin/Grafico_CurvasObs_Vs_CurvasSua_medellin_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - Medellin - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - Medellin - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_medellin/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - Medellin",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - Medellin",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_medellin/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",32),1:32)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8", "PO3", "PO4", "Pz", "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",32),1:32))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")


########################### EXEMPLARY PLOTS FDA ##########################
w <- "pre_alpha"
i <- 22
#for(i in 1:length(canales) ){

data_archivo <- archivo %>%
  dplyr::filter(channel == canales[23])%>%
  dplyr::select(subject, group, epoch, pre_alpha )



#Datos funcionales definitivos
data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = "pre_alpha" )
data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]



# par(mfrow = c(1,2) ) 
# matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
#         col = 1,
#         main = "HC",
#         xlab = "Epoch",
#         ylab = indicador[w] ,
#         xaxt = "n",
#         lty  = 1)
# 
# matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
#         col = 1,
#         main = "PD",
#         xlab = "Epoch",
#         ylab = indicador[w] ,
#         xaxt = "n",
#         lty  = 1)




data_wide_hc <- data_wide[data_wide$group=="hc",]
data_wide_hc <- data_wide_hc[, -c(1,2)]
data_wide_pd <- data_wide[data_wide$group=="pd",]
data_wide_pd <- data_wide_pd[, -c(1,2)]


TiempoEval = seq(1, dim(data_wide_hc)[2]) 
rangval    = c(1, dim(data_wide_hc)[2])

#### PARA MEDELLIN EL NÚMERO ÓPTIMO DE BASES ELEGIDO ES 21

fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 21, 4)


#fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
#fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)

fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )



#Extraer los datos funcionales
fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)





par(mfrow = c(2,3), mar = c(4, 4, 2,1))
matplot(t(data_wide_hc), type = "l", pch = 1,
        col = 1,
        main = "Observed curves - Non-PD - Medellin",
        xlab = "Epochs (5-sec)",
        ylab = "Pre-alpha relative PSD",
        ylim = c(0,0.7),
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)
matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
        col = 2,
        ylim = c(0,0.7),
        main = "Smoothed curves - Non-PD - Medellin",
        xlab = "Epochs (5-sec)",
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)


#Promedio funcional
Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)





plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
     main = "Functional mean - Non-PD - Medellin",
     ylab="",
     xlab = "Epochs (5-sec)",
     ylim = c(0,0.7),
     xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
lines(Prom_fd_HC_California, col = 2, lwd = 2)





matplot(t(data_wide_pd), type = "l", pch = 1,
        col = 1,
        main = "Observed curves - PD - Medellin",
        xlab = "Epochs (5-sec)",
        ylab = "Pre-alpha relative PSD",
        ylim = c(0,0.7),
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)

matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
        col = 2,
        ylim = c(0,0.7),
        main = "Smoothed curves - PD - Medellin",
        xlab = "Epochs (5-sec)",
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)




plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
     main = "Functional mean - PD - Medellin",
     xlab = "Epochs (5-sec)",
     ylab="",
     ylim = c(0,0.7),
     xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
lines(Prom_fd_PD_California, col = 2, lwd = 2)






#Curva t - funcional

### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                   fdata_data_wide_pd_adf$fd,
                                   argvals = 1:max(TiempoEval),
                                   plotres = F)


t_funcional_california$pvals.pts


par(mfrow = c(1,1))
plot(Prom_fd_HC_California, lty = 1, col = "blue",
     lwd = 2,
     main = "Functional Mean - Oz Channel",
     ylab = "Pre-alpha relative PSD",
     xlab = "Epochs (5-sec)",
     ylim = c(0,0.3),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
legend("bottomright", 
       c("Non-PD", "PD"),
       col = c("blue", "red"), lwd = c(2, 2))


###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "medellin_alpha_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_medellin_alpha_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/alpha.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nAlpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/medellin_beta_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_medellin_beta_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/beta.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nBeta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_higuchi_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/medellin_alpha_theta_channel_epoch.xlsx")
California_higuchi_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_medellin_alpha_theta_channel_epoch.xlsx")

California_higuchi_fd_p_long <-  California_higuchi_fd_p %>%
  mutate(Epoch = row.names(California_higuchi_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd_t_long <-  California_higuchi_fd_t %>%
  mutate(Epoch = row.names(California_higuchi_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd <- California_higuchi_fd_t_long %>%
  inner_join(y = California_higuchi_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/alpha_theta.png", width = 7, height = 5, units = "in", res = 300)
California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nAlpha/theta ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()





California_hjort_complexity_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/medellin_delta_channel_epoch.xlsx")
California_hjort_complexity_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_medellin_delta_channel_epoch.xlsx")

California_hjort_complexity_p_long <-  California_hjort_complexity_p %>%
  mutate(Epoch = row.names(California_hjort_complexity_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity_t_long <-  California_hjort_complexity_t %>%
  mutate(Epoch = row.names(California_hjort_complexity_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity <- California_hjort_complexity_t_long %>%
  inner_join(y = California_hjort_complexity_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/delta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nDelta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_hjort_mobility_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/medellin_theta_channel_epoch.xlsx")
California_hjort_mobility_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_medellin_theta_channel_epoch.xlsx")

California_hjort_mobility_p_long <-  California_hjort_mobility_p %>%
  mutate(Epoch = row.names(California_hjort_mobility_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility_t_long <-  California_hjort_mobility_t %>%
  mutate(Epoch = row.names(California_hjort_mobility_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility <- California_hjort_mobility_t_long %>%
  inner_join(y = California_hjort_mobility_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/theta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nTheta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()







California_katz_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/medellin_pre_alpha_channel_epoch.xlsx")
California_katz_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_medellin_pre_alpha_channel_epoch.xlsx")

California_katz_fd_p_long <-  California_katz_fd_p %>%
  mutate(Epoch = row.names(California_katz_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd_t_long <-  California_katz_fd_t %>%
  mutate(Epoch = row.names(California_katz_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd <- California_katz_fd_t_long %>%
  inner_join(y = California_katz_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/prealpha.png", width = 7, height = 5, units = "in", res = 300)
California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nPre-alpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()






fig1 = California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Delta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)




fig2 = California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Theta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)



fig3 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig4 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Beta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_bands_medellin   <- ggarrange(fig1, fig2,fig3, fig4, common.legend = T,
                                     legend = "right",
                                     ncol = 4, nrow = 1, 
                                     font.label = list(face = "bold", family = "mono"))

figure_bands_medellin   <- annotate_figure(figure_bands_medellin, top = text_grob("Medellin dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/bands_medellin.png", width = 35, height = 5, units = "in", res = 300)
figure_bands_medellin
dev.off()

fig5 = California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Pre-alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)

fig6 = California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha/theta relative PSD ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)















#####################   DOMINANT MEDELLIN  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "medellin"
archivo <- "df_medellin_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_med.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  setwd("D:/escritorio/tSNE_entropia_R/spectral/")
  
  pdf(paste0("ADF_medellin/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador
#######################################################################################
###############################################################################################################################
################################# FUNCTIONAL DATA ANALYSIS ################################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 20, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 20, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_medellin/Grafico_CurvasObs_Vs_CurvasSua_medellin_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - medellin",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - medellin",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_medellin/Grafico_CurvasObs_Vs_CurvasSua_medellin_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - medellin - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - medellin - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_medellin/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - medellin",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - medellin",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    max(TiempoEval)
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_medellin/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",32),1:32)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8", "PO3", "PO4", "Pz", "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",32),1:32))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")


















##################### PLOTS FDA #############
###############################################################################################  MEDELLIN  ##########################
###############################################################################



#Directorio de trabajo
#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/spectral")


carpeta <- "medellin"
archivo <- "nonorm_spectral_antioquia_5s_no_overlapping.csv"
indicador <- c("delta","theta","alpha"      
               ,"beta","pre_alpha","alpha_theta")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


# metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_med.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}


########################### EXEMPLARY PLOTS FDA ##########################
w <- "pre_alpha"
i <- 22
#for(i in 1:length(canales) ){

data_archivo <- archivo %>%
  dplyr::filter(channel == canales[22])%>%
  dplyr::select(subject, group, epoch, pre_alpha )



#Datos funcionales definitivos
data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = "pre_alpha" )
data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]



# par(mfrow = c(1,2) ) 
# matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
#         col = 1,
#         main = "HC",
#         xlab = "Epoch",
#         ylab = indicador[w] ,
#         xaxt = "n",
#         lty  = 1)
# 
# matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
#         col = 1,
#         main = "PD",
#         xlab = "Epoch",
#         ylab = indicador[w] ,
#         xaxt = "n",
#         lty  = 1)




data_wide_hc <- data_wide[data_wide$group=="hc",]
data_wide_hc <- data_wide_hc[, -c(1,2)]
data_wide_pd <- data_wide[data_wide$group=="pd",]
data_wide_pd <- data_wide_pd[, -c(1,2)]


TiempoEval = seq(1, dim(data_wide_hc)[2]) 
rangval    = c(1, dim(data_wide_hc)[2])

#### PARA MEDELLIN EL NÚMERO ÓPTIMO DE BASES ELEGIDO ES 21

fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 21, 4)


#fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
#fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)

fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )



#Extraer los datos funcionales
fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)



#Promedio funcional
Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)





#Curva t - funcional

### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                   fdata_data_wide_pd_adf$fd,
                                   argvals = 1:max(TiempoEval),
                                   plotres = F)


t_funcional_california$pvals.pts
#par(mfrow = c(2,2), mar = c(4, 4, 2,1))

###### GRID 1
plot(Prom_fd_HC_California, lty = 1, col = "blue",
     lwd = 2,
     main = "Functional Mean Medellin - Channel O2",
     ylab = "Pre-alpha relative PSD",
     xlab = "Epochs (5-sec)",
     ylim = c(0,0.3),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
legend("bottomright", 
       c("Non-PD", "PD"),
       col = c("blue", "red"), lwd = c(2, 2))

###### GRID 2
plot(Prom_fd_HC_California, lty = 1, col = "blue",
     lwd = 2,
     main = "Functional Mean Iowa - Channel O2",
     ylab = "Pre-alpha relative PSD",
     xlab = "Epochs (5-sec)",
     ylim = c(0,0.3),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)

###### GRID 3

plot(Prom_fd_HC_California, lty = 1, col = "blue",
     lwd = 2,
     main = "Functional Mean Finland - Channel O2",
     ylab = "Pre-alpha relative PSD",
     xlab = "Epochs (5-sec)",
     ylim = c(0,0.3),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)

###### GRID 4
plot(Prom_fd_HC_California, lty = 1, col = "blue",
     lwd = 2,
     main = "Functional Mean California - Channel O2 ",
     ylab = "Pre-alpha relative PSD",
     xlab = "Epochs (5-sec)",
     ylim = c(0,0.3),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)







#####################   DOMINANT MEDELLIN  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "medellin"
archivo <- "df_medellin_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_med.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  setwd("D:/escritorio/tSNE_entropia_R/spectral/")
  
  pdf(paste0("ADF_medellin/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador
#####################   DOMINANT IOWA  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "iowa"
archivo <- "df_iowa_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)

#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_iowa.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <- archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" ,"P4",  "P8" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))



#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador

#####################   DOMINANT FINLAND  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "finland"
archivo <- "df_finland_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_fin.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador

#####################   DOMINANT CALIFORNIA  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "california"
archivo <- "df_california_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_cal.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")
Iow_Par$subject <- as.numeric(Iow_Par$subject)

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador


########################### EXEMPLARY PLOTS FDA ##########################
w <- "df_ep"
i <- 21
#for(i in 1:length(canales) ){

data_archivo <- archivo %>%
  dplyr::filter(channel == canales[25])%>%
  dplyr::select(subject, group, epoch, df_ep )



#Datos funcionales definitivos
data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = "df_ep" )
data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]



# par(mfrow = c(1,2) ) 
# matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
#         col = 1,
#         main = "HC",
#         xlab = "Epoch",
#         ylab = indicador[w] ,
#         xaxt = "n",
#         lty  = 1)
# 
# matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
#         col = 1,
#         main = "PD",
#         xlab = "Epoch",
#         ylab = indicador[w] ,
#         xaxt = "n",
#         lty  = 1)




data_wide_hc <- data_wide[data_wide$group=="hc",]
data_wide_hc <- data_wide_hc[, -c(1,2)]
data_wide_pd <- data_wide[data_wide$group=="pd",]
data_wide_pd <- data_wide_pd[, -c(1,2)]


TiempoEval = seq(1, dim(data_wide_hc)[2]) 
rangval    = c(1, dim(data_wide_hc)[2])

#### PARA MEDELLIN EL NÚMERO ÓPTIMO DE BASES ELEGIDO ES 21

fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 21, 4)


#fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
#fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)

fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )



#Extraer los datos funcionales
fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)





par(mfrow = c(2,3), mar = c(4, 4, 2,1))
matplot(t(data_wide_hc), type = "l", pch = 1,
        col = 1,
        main = "Observed curves - Non-PD - Medellin",
        xlab = "Epochs (5-sec)",
        ylab = "Dominant Frequency (Hz)",
        ylim = c(4, 16),
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)
matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
        col = 2,
        ylim = c(4, 16),
        main = "Smoothed curves - Non-PD - Medellin",
        xlab = "Epochs (5-sec)",
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)


#Promedio funcional
Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)





plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
     main = "Functional mean - Non-PD - Medellin",
     ylab="",
     xlab = "Epochs (5-sec)",
     ylim = c(4, 16),
     xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
lines(Prom_fd_HC_California, col = 2, lwd = 2)





matplot(t(data_wide_pd), type = "l", pch = 1,
        col = 1,
        main = "Observed curves - PD - Medellin",
        xlab = "Epochs (5-sec)",
        ylab = "Dominant Frequency (Hz)",
        ylim = c(4, 16),
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)

matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
        col = 2,
        ylim = c(4, 16),
        main = "Smoothed curves - PD - Medellin",
        xlab = "Epochs (5-sec)",
        xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
        lty = 1,
        lwd = 2)




plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
     main = "Functional mean - PD - Medellin",
     xlab = "Epochs (5-sec)",
     ylab="",
     ylim = c(4, 16),
     xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
lines(Prom_fd_PD_California, col = 2, lwd = 2)






#Curva t - funcional

### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                   fdata_data_wide_pd_adf$fd,
                                   argvals = 1:max(TiempoEval),
                                   plotres = F)


t_funcional_california$pvals.pts

###### PLOTS GRID 1

par(mfrow = c(1,1))
fig_med = plot(Prom_fd_HC_California, lty = 1, col = "blue",
               lwd = 2,
               main = "Functional Mean Medellin - P4 Channel",
               ylab = "Dominant Frequency (Hz)",
               xlab = "Epochs (5-sec)",
               ylim = c(7,11),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)

##### PLOTS GRID 2

par(mfrow = c(1,1))
fig_iow = plot(Prom_fd_HC_California, lty = 1, col = "blue",
               lwd = 2,
               main = "Functional Mean Iowa - P4 Channel",
               ylab = "Dominant Frequency (Hz)",
               xlab = "Epochs (5-sec)",
               ylim = c(7,11),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)


##### PLOTS GRID 3

par(mfrow = c(1,1))
fig_fin = plot(Prom_fd_HC_California, lty = 1, col = "blue",
               lwd = 2,
               main = "Functional Mean Finland - P4 Channel",
               ylab = "Dominant Frequency (Hz)",
               xlab = "Epochs (5-sec)",
               ylim = c(7,11),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)

##### PLOTS GRID 4

par(mfrow = c(1,1))
fig_cal = plot(Prom_fd_HC_California, lty = 1, col = "blue",
               lwd = 2,
               main = "Functional Mean California - P4 Channel",
               ylab = "Dominant Frequency (Hz)",
               xlab = "Epochs (5-sec)",
               ylim = c(7,11),)
lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)



par(mfrow = c(2,3), mar = c(4, 4, 2,1))
fig_med
fig_iow
fig_fin
fig_cal


legend("bottomright", 
       c("Non-PD", "PD"),
       col = c("blue", "red"), lwd = c(2, 2))



###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "medellin_dfv_ch_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_medellin_dfv_ch_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/dfv.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nDominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/medellin_df_ep_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_medellin_df_ep_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/df.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Medellin\nDominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




fig8 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig7 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_spectral_medellin   <- ggarrange(fig6, fig5,fig7, fig8, common.legend = T,
                                        legend = "right",
                                        ncol = 4, nrow = 1, 
                                        font.label = list(face = "bold", family = "mono"))

figure_spectral_medellin   <- annotate_figure(figure_spectral_medellin, top = text_grob("Medellin dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_medellin")
png(filename = "./HEATMAPS/spectral.png", width = 35, height = 5, units = "in", res = 300)
figure_spectral_medellin
dev.off()


































############################################################################################  CALIFORNIA  ##########################
###############################################################################



#Directorio de trabajo
#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/spectral")


carpeta <- "california"
archivo <- "nonorm_spectral_california_5s_no_overlapping.csv"
indicador <- c("delta","theta","alpha"      
               ,"beta","pre_alpha","alpha_theta")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)



archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  
  pdf(paste0("ADF_california/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}
#termina la corrida para cada indicador#######################################################################################
##### DATOS FUNCIONALES PARA CADA MÉTRICA ####################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo




for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    #### PARA MEDELLIN EL NÚMERO ÓPTIMO DE BASES ELEGIDO ES 30
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_california/Grafico_CurvasObs_Vs_CurvasSua_california_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - California",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - California",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_california/Grafico_CurvasObs_Vs_CurvasSua_california_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - California - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - California - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_california/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - California",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - California",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_california/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",32),1:32)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8", "PO3", "PO4", "Pz", "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",32),1:32))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")



###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "california_alpha_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_california_alpha_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/alpha.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nAlpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/california_beta_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_california_beta_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/beta.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nBeta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_higuchi_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/california_alpha_theta_channel_epoch.xlsx")
California_higuchi_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_california_alpha_theta_channel_epoch.xlsx")

California_higuchi_fd_p_long <-  California_higuchi_fd_p %>%
  mutate(Epoch = row.names(California_higuchi_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd_t_long <-  California_higuchi_fd_t %>%
  mutate(Epoch = row.names(California_higuchi_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_higuchi_fd <- California_higuchi_fd_t_long %>%
  inner_join(y = California_higuchi_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/alpha_theta.png", width = 7, height = 5, units = "in", res = 300)
California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nAlpha/theta ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()





California_hjort_complexity_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/california_delta_channel_epoch.xlsx")
California_hjort_complexity_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_california_delta_channel_epoch.xlsx")

California_hjort_complexity_p_long <-  California_hjort_complexity_p %>%
  mutate(Epoch = row.names(California_hjort_complexity_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity_t_long <-  California_hjort_complexity_t %>%
  mutate(Epoch = row.names(California_hjort_complexity_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_complexity <- California_hjort_complexity_t_long %>%
  inner_join(y = California_hjort_complexity_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/delta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nDelta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_hjort_mobility_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/california_theta_channel_epoch.xlsx")
California_hjort_mobility_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_california_theta_channel_epoch.xlsx")

California_hjort_mobility_p_long <-  California_hjort_mobility_p %>%
  mutate(Epoch = row.names(California_hjort_mobility_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility_t_long <-  California_hjort_mobility_t %>%
  mutate(Epoch = row.names(California_hjort_mobility_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_hjort_mobility <- California_hjort_mobility_t_long %>%
  inner_join(y = California_hjort_mobility_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/theta.png", width = 7, height = 5, units = "in", res = 300)
California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nTheta relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()







California_katz_fd_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/california_pre_alpha_channel_epoch.xlsx")
California_katz_fd_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_california_pre_alpha_channel_epoch.xlsx")

California_katz_fd_p_long <-  California_katz_fd_p %>%
  mutate(Epoch = row.names(California_katz_fd_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd_t_long <-  California_katz_fd_t %>%
  mutate(Epoch = row.names(California_katz_fd_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_katz_fd <- California_katz_fd_t_long %>%
  inner_join(y = California_katz_fd_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))

# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/prealpha.png", width = 7, height = 5, units = "in", res = 300)
California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nPre-alpha relative power spectral density") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()






fig1 = California_hjort_complexity %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Delta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)




fig2 = California_hjort_mobility %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Theta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel",legend.position="none", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)



fig3 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig4 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Beta relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


figure_bands_california   <- ggarrange(fig1, fig2,fig3, fig4, common.legend = T,
                                       legend = "right",
                                       ncol = 4, nrow = 1, 
                                       font.label = list(face = "bold", family = "mono"))

figure_bands_california   <- annotate_figure(figure_bands_california, top = text_grob("California dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/bands_california.png", width = 20, height = 5, units = "in", res = 300)
figure_bands_california
dev.off()

fig5 = California_katz_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Pre-alpha relative PSD") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)

fig6 = California_higuchi_fd %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Alpha/theta relative PSD ratio") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)















#####################   DOMINANT CALIFORNIA  ############################



#Directorio de trabajo
setwd("D:/escritorio/tSNE_entropia_R/dominant/")


carpeta <- "california"
archivo <- "df_california_5s_no_overlapping.csv"
indicador <- c("df_ep","dfv_ch")



archivo <- fread(paste0(carpeta,"/",archivo),
                 sep = ",",
                 dec = ".",
                 header = TRUE)


#Iowa metadata
Iow_Par <- read.table(file = "D:/escritorio/tSNE_entropia_R/participants_cal.tsv", sep = "\t", header = TRUE)
Iow_Par <- Iow_Par %>%
  separate(participant_id, c("Sub", "subject"), "-")

Iow_Par <- Iow_Par %>%
  dplyr::select(subject, group, age, gender)
#  select(subject, group, subgroup, moca_total, age)

archivo <- archivo %>%
  inner_join(Iow_Par, by="subject")

archivo <-archivo %>% filter_at(vars(channel), any_vars(. %in%  c("Fp1", "Fp2","AF3","AF4","F7","F3","Fz","F4", "F8" ,"FC5", "FC1", "FC2", "FC6", "T7",  "C3",  "Cz",  "C4",  "T8",  "CP5", "CP1", "CP2" ,"CP6" ,"P7" , "P3" , "Pz"  ,"P4",  "P8" , "PO3" ,"PO4" ,"O1" , "Oz" , "O2" )))                                                              


canales <- names(table(archivo$channel))


#w <- 1

for(w in 1:length(indicador)){#empieza la corrida para cada k
  
  #Cargar datos
  
  #archivo %>%
  #  ggplot(aes(x = epoch, y = permutation, color = subject)) + 
  #  geom_point()+
  #  geom_line()+
  #  facet_wrap( ~ group)
  
  
  data_archivo <- archivo %>%
    dplyr::filter(channel == canales[1])%>%
    dplyr::select(subject, group, epoch, indicador[w] )
  
  
  
  #data_wide <- impute(data_prueba$permutation, median)
  data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
  data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
  
  
  
  # par(mfrow = c(1,2) ) 
  # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "HC",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  # 
  # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
  #         col = 1,
  #         main = "PD",
  #         xlab = "Epoch",
  #         ylab = indicador[w] ,
  #         xaxt = "n",
  #         lty  = 1)
  
  
  
  
  data_wide_hc <- data_wide[data_wide$group=="hc",]
  data_wide_hc <- data_wide_hc[, -c(1,2)]
  data_wide_pd <- data_wide[data_wide$group=="pd",]
  data_wide_pd <- data_wide_pd[, -c(1,2)]
  
  #Numero de bases
  nbases <- seq(4,dim(data_wide_hc)[2],1)
  #lam    <- seq(0.1,1, 0.1)
  
  TiempoEval = seq(1, dim(data_wide_hc)[2]) 
  rangval    = c(1, dim(data_wide_hc)[2])
  
  #Seleccion del numero de bases
  
  #fdata_data_wide_hc <- fdata(data_wide_hc, TiempoEval, rangval)
  #k_hc <- optim.basis(fdata_data_wide_hc, lambda = lam, numbasis = nbases)
  #k_hc$numbasis.opt
  
  #fdata_data_wide_pd <- fdata(data_wide_pd, TiempoEval, rangval)
  #k_pd <- optim.basis(fdata_data_wide_pd, lambda = lam, numbasis = nbases)
  #k_pd$numbasis.opt
  
  k <- 5:(dim(data_wide_hc)[2]-2)
  #knots    = c(seq(0,dim(data_wide_hc)[2],1))
  
  Sesgo_hc <- NULL
  Sesgo_pd <- NULL
  Sesgo_Prom_hc <- NULL
  Sesgo_Prom_pd <- NULL
  
  #Sesgo_Prom_hc <- NULL
  #Sesgo_Prom_pd <- NULL
  
  
  
  
  
  for(j in 1:length(k)){
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = k[j], 4)
    
    
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    #RECM_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)^2
    #RECM_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)^2
    
    #Sesgo_hc <- (t(fvalores_data_wide_hc_adf) - data_wide_hc)
    #Sesgo_pd <- (t(fvalores_data_wide_pd_adf) - data_wide_pd)
    
    for(i in 1:dim(data_wide_hc)[1]){Sesgo_hc[i] <-  bias( as.numeric(t(fvalores_data_wide_hc_adf)[i,]), as.numeric(data_wide_hc[i,]) )}
    for(i in 1:dim(data_wide_pd)[1]){Sesgo_pd[i] <-  bias( as.numeric(t(fvalores_data_wide_pd_adf)[i,]), as.numeric(data_wide_pd[i,]) )}
    
    
    
    
    #matplot(RECM_hc, type = "l", col = "grey", ylim = c(-0.02, 0.02))
    #abline(h = 0, col = 2, lty = 2)
    
    Sesgo_Prom_hc[j] <-  mean(Sesgo_hc) 
    Sesgo_Prom_pd[j] <-  mean(Sesgo_pd) 
    
    #Sesgo_Prom_hc[i] <-  mean(apply(Sesgo_hc , 1, mean)) 
    #Sesgo_Prom_pd[i] <-  mean(apply(Sesgo_pd , 1, mean)) 
    
    
  }
  setwd("D:/escritorio/tSNE_entropia_R/spectral/")
  
  pdf(paste0("ADF_california/Sesgo_promedio",indicador[w], ".pdf"))
  plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
       main = indicador[w],
       xlab = "Basis number (k)",
       ylab = "Mean Bias",
       ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd) )) 
  lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  abline(h = 0, col = "red", lty = 2)
  legend("bottomright", 
         c("HC", "PD"),
         col = c(1, 2), lwd = c(2, 2))
  dev.off()
  
  
  # pdf(paste0("ADF_California/Sesgo_promedio_",indicador[w], ".pdf"))
  # plot(k,  Sesgo_Prom_hc, type = "b", pch = 19,
  #      main = indicador[w],
  #      xlab = "Num bases (k)",
  #      ylab = "Sesgo promedio",
  #      ylim = c(min(Sesgo_Prom_hc, Sesgo_Prom_pd), max(Sesgo_Prom_hc, Sesgo_Prom_pd )))
  # lines(k, Sesgo_Prom_pd, type = "b", pch = 19, col = 2)
  # abline(h = 0, col = "red", lty = 2)
  # dev.off()
  
  #plot(k, Sesgo_Prom_hc)
  #lines(k, RECM_Prom_hc)
  
  #plot(Sesgo_Prom_hc,RECM_Prom_hc, type = "l")
  #abline(v = 0, col = "red", lty = 2)
  
  
}#termina la corrida para cada indicador
#######################################################################################
###############################################################################################################################
################################# FUNCTIONAL DATA ANALYSIS ################################################################


per_dif       <- NULL
p_value_epoch <- list()
aux <- list()
t_value_epoch <- list()#nuevo
aux_t_values <- list()#nuevo





for(w in 1:length(indicador) ){
  for(i in 1:length(canales) ){
    
    data_archivo <- archivo %>%
      dplyr::filter(channel == canales[i])%>%
      dplyr::select(subject, group, epoch, indicador[w] )
    
    
    
    #Datos funcionales definitivos
    data_wide <- dcast(data_archivo, subject + group ~ epoch, value.var = indicador[w] )
    data_wide <- data_wide[,1:length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])]
    
    
    
    # par(mfrow = c(1,2) ) 
    # matplot(t(data_wide[data_wide$group=="hc",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "HC",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    # 
    # matplot(t(data_wide[data_wide$group=="pd",]), type = "l", pch = 1,
    #         col = 1,
    #         main = "PD",
    #         xlab = "Epoch",
    #         ylab = indicador[w] ,
    #         xaxt = "n",
    #         lty  = 1)
    
    
    
    
    data_wide_hc <- data_wide[data_wide$group=="hc",]
    data_wide_hc <- data_wide_hc[, -c(1,2)]
    data_wide_pd <- data_wide[data_wide$group=="pd",]
    data_wide_pd <- data_wide_pd[, -c(1,2)]
    
    
    TiempoEval = seq(1, dim(data_wide_hc)[2]) 
    rangval    = c(1, dim(data_wide_hc)[2])
    
    
    
    fdata_data_wide_hc_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
    fdata_data_wide_pd_bases <- create.bspline.basis(rangval, nbasis = 21, 4)
    
    
    #fdata_data_wide_hc_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    #fdata_data_wide_pd_bases  = create.bspline.basis(rangval, nbasis = 15, norder = 4)
    
    fdata_data_wide_hc_adf = smooth.basis(TiempoEval,t(data_wide_hc), fdata_data_wide_hc_bases )
    fdata_data_wide_pd_adf = smooth.basis(TiempoEval,t(data_wide_pd), fdata_data_wide_pd_bases )
    
    
    
    #Extraer los datos funcionales
    fvalores_data_wide_hc_adf = eval.fd(TiempoEval, fdata_data_wide_hc_adf$fd)
    fvalores_data_wide_pd_adf = eval.fd(TiempoEval, fdata_data_wide_pd_adf$fd)
    
    
    
    
    pdf(paste0("ADF_california/Grafico_CurvasObs_Vs_CurvasSua_california_HC_", indicador[w], ".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_hc), type = "l", pch = 1,
            col = 1,
            main = "Curvas observadas - HC - california",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_hc_adf, type = "l", pch = 1,
            col = 2,
            main = "Curvas suavizadas - HC - california",
            xlab = "Epochs (5-sec)",
            ylab = indicador[w],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    
    
    pdf(paste0("ADF_california/Grafico_CurvasObs_Vs_CurvasSua_california_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2), mar = c(4, 4, 2,2))
    matplot(t(data_wide_pd), type = "l", pch = 1,
            col = 1,
            main = paste0("Curvas observadas - PD - california - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    
    matplot(fvalores_data_wide_pd_adf, type = "l", pch = 1,
            col = 2,
            main = paste0("Curvas suavizadas - PD - california - ", indicador[w]),
            xlab = "Epochs (5-sec)",
            ylab = indicador[1],
            xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])),
            lty = 1,
            lwd = 2)
    dev.off()
    
    
    
    #Promedio funcional
    Prom_fd_HC_California <- mean.fd(fdata_data_wide_hc_adf$fd)
    DE_fd_HC_California   <- std.fd(fdata_data_wide_hc_adf$fd)
    Prom_fd_PD_California <- mean.fd(fdata_data_wide_pd_adf$fd)
    DE_fd_PD_California   <- std.fd(fdata_data_wide_pd_adf$fd)
    
    
    
    pdf(paste0("ADF_california/DatosFuncionales_PromFuncional_HC_PD_",indicador[w],".pdf"))
    par(mfrow = c(1,2))
    plot(fdata_data_wide_hc_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - HC - california",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])) )
    lines(Prom_fd_HC_California, col = 2, lwd = 2)
    
    
    plot(fdata_data_wide_pd_adf, lty = 1, col = "grey60",
         main = "Curvas suavizadas - PD - california",
         xlab = "Epochs (5-sec)",
         ylab = indicador[w],
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),
         xlim = c(0, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0])))
    lines(Prom_fd_PD_California, col = 2, lwd = 2)
    
    dev.off()
    
    
    
    max(TiempoEval)
    
    #Curva t - funcional
    
    ### NÚMERO MÁXIMO DE ÉPOCAS ES EL ARGVALS MAYOR
    t_funcional_california <- tperm.fd(fdata_data_wide_hc_adf$fd, 
                                       fdata_data_wide_pd_adf$fd,
                                       argvals = 1:max(TiempoEval),
                                       plotres = F)
    
    
    t_funcional_california$pvals.pts
    
    
    par(mfrow = c(1,2))
    plot(Prom_fd_HC_California, lty = 1, col = "blue",
         lwd = 2,
         main = "Functional mean",
         ylab = indicador[w],
         xlab = "Epochs (5-sec)",
         ylim = c(min(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y), max(fdata_data_wide_hc_adf$y, fdata_data_wide_pd_adf$y)),)
    lines(Prom_fd_PD_California, lty = 1, col = "red", lwd = 2)
    legend("bottomright", 
           c("HC", "PD"),
           col = c("blue", "red"), lwd = c(2, 2))
    
    
    pdf(paste0("ADF_california/DatosFuncionales_Valorp_ttest_",indicador[w],"_Canal_",i,".pdf"))
    plot(seq(1, length(colSums(is.na(data_wide[,-c(1,2)]))[colSums(is.na(data_wide[,-c(1,2)])) == 0]), length=length(t_funcional_california $pvals.pts)),
         t_funcional_california$pvals.pts,   xaxt = "n", type = "b",
         pch = 19,
         ylim = c(0,1),
         xlab = "Epoch",
         ylab = "p value")
    abline(h = 0.05, col = 2, lty = 2)
    dev.off()
    
    p_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$pvals.pts)
    t_value_epoch[[i]] <- cbind.data.frame(t_funcional_california$Tvalores) #nuevo
    
    #per_dif[i] <-  length(t_funcional_california$pvals.pts[t_funcional_california$pvals.pts<0.05])/length(t_funcional_california$pvals.pts)
  }
  
  
  aux[[w]] <- do.call(cbind.data.frame, p_value_epoch)
  colnames(aux[[w]]) <- paste0(rep("Channel_",32),1:32)
  write.xlsx(aux[[w]], paste0(carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )
  
  
  aux_t_values[[w]] <- do.call(cbind.data.frame, t_value_epoch)#nuevo
  colnames(aux_t_values[[w]]) <- paste0(rep("Channel_",length(canales)),1:length(canales))#nuevo
  write.xlsx(aux_t_values[[w]], paste0("tvalue_", carpeta[1],"_",indicador[w],"_channel_epoch.xlsx"), overwrite = T )#nuevo
  
}




nombre.canales <- cbind.data.frame("Nombre" =  c("AF3", "AF4", "C3",  "C4",  "CP1", "CP2", "CP5", "CP6", "Cz",  "F3",  "F4",  "F7",  "F8",
                                                 "FC1", "FC2", "FC5", "FC6", "Fp1", "Fp2", "Fz",  "O1", "O2",  "Oz",  "P3",  "P4",
                                                 "P7",  "P8", "PO3", "PO4", "Pz", "T7",  "T8"),
                                   "id" = paste0(rep("Channel_",32),1:32))

library(openxlsx)
write.xlsx(nombre.canales, "nombre.canales.xlsx")










###### HEATMAPS ####################


setwd("D:/escritorio/tSNE_entropia_R/spectral")

# Lectura de CANALES

Canales <- read.table(file = "D:/escritorio/tSNE_entropia_R/spectral/Canales.csv", header = TRUE, sep = ",")
Canales <- within(data = Canales, expr = {
  Name <- factor(Name, levels = rev(Name), labels = rev(Name))
})
Canales$Name


California_approximate_p      <- read_excel(path = "california_dfv_ch_channel_epoch.xlsx")
California_approximate_t     <- read_excel(path = "tvalue_california_dfv_ch_channel_epoch.xlsx")

California_approximate_p_long <-  California_approximate_p %>%
  mutate(Epoch = row.names(California_approximate_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate_t_long <-  California_approximate_t %>%
  mutate(Epoch = row.names(California_approximate_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_approximate <- California_approximate_t_long %>%
  inner_join(y = California_approximate_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/dfv.png", width = 7, height = 5, units = "in", res = 300)
California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nDominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




California_detrended_fluct_p      <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/california_df_ep_channel_epoch.xlsx")
California_detrended_fluct_t     <- read_excel(path = "D:/escritorio/tSNE_entropia_R/spectral/tvalue_california_df_ep_channel_epoch.xlsx")

California_detrended_fluct_p_long <-  California_detrended_fluct_p %>%
  mutate(Epoch = row.names(California_detrended_fluct_p)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "P_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct_t_long <-  California_detrended_fluct_t %>%
  mutate(Epoch = row.names(California_detrended_fluct_t)) %>%
  pivot_longer(cols = starts_with("Channel_"), names_to = "Channel", names_prefix = "Channel_", values_to = "t_values") %>%
  mutate(Epoch = parse_number(Epoch), Channel = parse_number(Channel)) %>%
  inner_join(y = Canales, by = "Channel")

California_detrended_fluct <- California_detrended_fluct_t_long %>%
  inner_join(y = California_detrended_fluct_p_long, by = c("Channel", "Index", "Region", "Name", "Epoch")) %>%
  mutate(star = ifelse(P_values < 0.05, "*", ""))


# GrÃ¡ficos
setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/df.png", width = 7, height = 5, units = "in", res = 300)
California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "California\nDominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(panel.grid.minor.y = element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)
dev.off()




fig8 = California_approximate %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency Variability (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)),
        plot.caption.position = "panel", legend.position="none",  panel.grid.minor.y = element_blank(), strip.text.y=element_blank(),
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)


fig7 =California_detrended_fluct %>%
  mutate(star = ifelse(P_values < 0.05, "*", "")) %>%
  ggplot(aes(x = as.factor(Epoch), y = Name, fill = t_values)) +
  labs(y = "", x = "Epoch")  +
  geom_tile() + labs(caption = "Dominant Frequency (Hz)") +
  guides(fill = guide_colourbar(barwidth = 0.5,
                                barheight = 20))+
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))+   facet_grid(Region ~ ., space = 'free_y', scales = 'free_y', switch = 'y') +
  theme_classic(base_family = 'mono') +
  theme(plot.caption = element_text(vjust = 0, hjust = 0.5, size = 14, margin = margin(10, 0, 0, 0)), 
        plot.caption.position = "panel", panel.grid.minor.y = element_blank(), strip.text.y=element_blank(), 
        panel.spacing.y = unit(0,"line"), 
        strip.placement = 'outside', strip.background.y = element_blank()) + 
  scale_fill_gradientn(name = "t value", colors = viridis::cividis(n = 256, begin = 0, end=1),  limits=c(-5.73,5.73)) +
  #  scale_fill_gradientn(name = "t value", colors = rainbow(n = 10, start = 0, end = 0.7), limits=c(0,1)) +
  geom_text(aes(label = star), color = "black", size = 5, nudge_y = -0.2)

figure_spectral_california   <- ggarrange(fig6, fig5,fig7, fig8, common.legend = T,
                                          legend = "right",
                                          ncol = 4, nrow = 1, 
                                          font.label = list(face = "bold", family = "mono"))

figure_spectral_california   <- annotate_figure(figure_spectral_california, top = text_grob("California dataset", face = "bold", family = 'mono', size = 14))

setwd("D:/escritorio/tSNE_entropia_R/spectral/ADF_california")
png(filename = "./HEATMAPS/spectral.png", width = 20, height = 5, units = "in", res = 300)
figure_spectral_california
dev.off()
















































################################################################################################################## PAPER FIGURES #############


figure_spectral_all   <- ggarrange(figure_spectral_california, figure_spectral_finland,figure_spectral_iowa, figure_spectral_medellin,
                                   ncol = 1, nrow = 4, common.legend = T, legend = "right",
                                   font.label = list(face = "bold", family = "mono"))

setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "spectral_all.png", width = 35, height = 22, units = "in", res = 300)
figure_spectral_all
dev.off()




figure_bands_all   <- ggarrange(figure_bands_california, figure_bands_finland,figure_bands_iowa, figure_bands_medellin,
                                ncol = 1, nrow = 4, common.legend = T, legend = "right",
                                font.label = list(face = "bold", family = "mono"))

setwd("D:/escritorio/tSNE_entropia_R/Publication")
png(filename = "bands_all.png", width = 35, height = 22, units = "in", res = 300)
figure_bands_all
dev.off()