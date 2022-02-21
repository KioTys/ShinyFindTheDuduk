Registres <- function(f, freqs, labels) {
  freqs  <- c(0, freqs, 1e5)
  labels <- c(labels, "Out_of_range")
  class  <- 0
  i      <- 1
  while(class == 0 & i<length(freqs)) {
    if (f>=freqs[i] & f<freqs[i+1]) {
      class <- 1
      return (labels[i])
      i <- i+1
    }
    else {
      i <- i+1
    }
  }
}




PredErr <- function(X, Y) {
  mean(X!=Y)
}


PrepMorceau <- function(morceau, InstrInput) {
  library(tidyverse)
  library(seewave)
  library(tuneR)
  
  sound   <- readMP3(morceau)
  
  freqs  <- c(60, 200, 640, 1625, 7000, 10000) # limites: 0-60, registre "Extremement grave", 60-200, registre "Grave", etc.
  labels <- c("ExtrGrave", "Grave", "BasMed", "HautsMed", "Aigu", "ExtrAigu")
  n.labels <- length(labels)
  
  data_final <- as.data.frame(matrix(ncol = 1+1+1+1+3*length(labels)+1, nrow=1)) # On cree un df avecune ligne de NA, qu'il faudra supprimer a la fin
  names(data_final) <- c("Morceau",
                         "t0",
                         "Duree_initiale",
                         "Duree_finale",
                         paste("A_moy", labels, sep="_"),
                         paste("A_q05", labels, sep="_"),
                         paste("A_q95", labels, sep="_"),
                         "Instrument")
  
  f         <- sound@samp.rate
  d1        <- length(sound@left)/f 
  
  nt    <- d1*f                # Nombre de points dans l'echantillon temporel
  nf    <- nt/2                # Nombre de points dans l'echantillon fr?quentiel
  
  ti    <- 20                  # t0 pour le decoupage des morceaux (20 secondes pour les applaudissements, etc.)
  
  data  <- data.frame(matrix(nrow = 1, ncol=1+1+1+1+3*length(labels)+1))  # Source, plage, d1, d2, meanS, Q05S, Q95S, Instrument
  
  sound_c_f   <- as.data.frame(seewave::spec(sound, f = f, plot = F, fftw = T))  # D?composition en frequences et stockage dans un df
  sound_c_f$Registre <- unlist(apply(data.frame(sound_c_f$x*1e3), FUN=Registres, MARGIN=1, freqs, labels))
  
  out <- sound_c_f %>% 
    group_by(Registre) %>% 
    summarise(mean=mean(y), q05=quantile(y, 0.05), q95=quantile(y, 0.95))
  
  data[1, 1]                                              <- morceau
  data[1, 2]                                              <- ti
  data[1, 3:4]                                            <- c(d1, d2)
  data[1, 5:(5+length(labels)-1)]                         <- out[1,-1]
  data[1, (5+length(labels)):(5+2*length(labels)-1)]      <- out[2,-1]
  data[1, (5+2*length(labels)):(5+3*length(labels)-1)]    <- out[3,-1]
  data[1, 5+3*length(labels)]                             <- InstrInput  
  
  data <- as.data.frame(data)
  names(data) <- c("Morceau",
                   "t0",
                   "Duree_initiale",
                   "Duree_finale",
                   paste("A_moy", labels, sep="_"),
                   paste("A_q05", labels, sep="_"),
                   paste("A_q95", labels, sep="_"),
                   "Instrument")  
  
  return(data)
}

PrepMorceau_Cut <- function(morceau, InstrInput, seq) {
  library(tidyverse)
  library(seewave)
  library(tuneR)
  
  sound   <- readMP3(morceau)
  spectro(sound, sound@samp.rate, 
          flim = c(0,6), 
          flog = T, 
          fastdisp = T, 
          osc = T, 
          grid = T,
          tlab = "Temps (s)",
          flab = "Frequence (kHz)")
  
  freqs  <- c(60, 200, 640, 1625, 7000, 10000) # limites: 0-60, registre "Extremement grave", 60-200, registre "Grave", etc.
  labels <- c("ExtrGrave", "Grave", "BasMed", "HautsMed", "Aigu", "ExtrAigu")
  n.labels <- length(labels)
  
  d_seq <- seq
  
  data_final <- as.data.frame(matrix(ncol = 1+1+1+1+3*length(labels)+1, nrow=1)) # On cr?e un df avecune ligne de NA, qu'il faudra supprimer ? la fin
  names(data_final) <- c("Morceau",
                         "t0",
                         "Duree_initiale",
                         "Duree_finale",
                         paste("A_moy", labels, sep="_"),
                         paste("A_q05", labels, sep="_"),
                         paste("A_q95", labels, sep="_"),
                         "Instrument")
  
  f         <- sound@samp.rate
  d1        <- length(sound@left)/f # Duree initiale
  # sound     <- zapsilw(sound, threshold = 1, output = "Wave", f = f,plot = F) # Suppression des silences - Semble assez chronophage
  d2        <- length(sound@left)/f # Duree finale (apr?s suppression des silences)
  
  nt    <- d_seq*f                # Nombre de points dans l'?chantillon temporel
  nf    <- nt/2                   # Nombre de points dans l'?chantillon fr?quentiel
  
  ti    <- d_seq                  # t0 pour le decoupage des morceaux (20 secondes pour les applaudissements, etc.)
  
  data  <- data.frame(matrix(nrow = (d2-ti)%/%d_seq, ncol=1+1+1+1+3*length(labels)+1))  # Source, plage, d1, d2, meanS, Q05S, Q95S, Instrument
  
  j     <- 1
  if(d2<1e9) { # Pour pouvoir eventuellement travailler uniquement sur les "petits" morceaux
    while ((ti+d_seq)<d2) { 
      print(paste(ti, ti + d_seq, d2, dim(data_final)[1]+j-1))    
      
      sound_c     <- cutw(sound, from = ti, to = ti + d_seq, f = f)                    # Extraction d'une s?quence de dur?e d_seq
      sound_c_f   <- as.data.frame(seewave::spec(sound_c, f = f, plot = F, fftw = T))  # D?composition en fr?quences et stockage dans un df
      sound_c_f$Registre <- unlist(apply(data.frame(sound_c_f$x*1e3), FUN=Registres, MARGIN=1, freqs, labels))
      
      if(max(abs(sound_c))>0) {
        out <- sound_c_f %>% 
          group_by(Registre) %>% 
          summarise(mean=mean(y), q05=quantile(y, 0.05), q95=quantile(y, 0.95))
        
        data[j, 1]                                              <- morceau
        data[j, 2]                                              <- ti
        data[j, 3:4]                                            <- c(d1, d2)
        data[j, 5:(5+length(labels)-1)]                         <- out[1,-1]
        data[j, (5+length(labels)):(5+2*length(labels)-1)]      <- out[2,-1]
        data[j, (5+2*length(labels)):(5+3*length(labels)-1)]    <- out[3,-1]
        data[j, 5+3*length(labels)]                             <- InstrInput  
        
        j <- j+1
        
      }
      
      ti <- ti + d_seq
      
      data <- as.data.frame(data)
      names(data) <- c("Morceau",
                       "t0",
                       "Duree_initiale",
                       "Duree_finale",
                       paste("A_moy", labels, sep="_"),
                       paste("A_q05", labels, sep="_"),
                       paste("A_q95", labels, sep="_"),
                       "Instrument")
    }
  }
  return(na.omit(rbind(data_final, data)))
}

PrepMorceau_Cut_Overlap <- function(morceau, InstrInput, d_seq, overlap=d_seq/2) {
  library(tidyverse)
  library(seewave)
  library(tuneR)
  
  sound   <- readMP3(morceau)
  spectro(sound, sound@samp.rate, 
          flim = c(0,6), 
          flog = T, 
          fastdisp = T, 
          osc = T, 
          grid = T,
          tlab = "Temps (s)",
          flab = "Frequence (kHz)")
  
  freqs  <- c(60, 200, 640, 1625, 7000, 10000) # limites: 0-60, registre "Extremement grave", 60-200, registre "Grave", etc.
  labels <- c("ExtrGrave", "Grave", "BasMed", "HautsMed", "Aigu", "ExtrAigu")
  n.labels <- length(labels)
  
  data_final <- as.data.frame(matrix(ncol = 1+1+1+1+1+3*length(labels)+1, nrow=1)) # On cree un df avec une ligne de NA, qu'il faudra supprimer ? la fin
  names(data_final) <- c("Morceau",
                         "t0",
                         "Duree_initiale",
                         "Duree_finale",
                         paste("A_moy", labels, sep="_"),
                         paste("A_q05", labels, sep="_"),
                         paste("A_q95", labels, sep="_"),
                         "Instrument",
                         "Temps")
  
  f         <- sound@samp.rate
  d1        <- length(sound@left)/f # Duree initiale
  # sound     <- zapsilw(sound, threshold = 1, output = "Wave", f = f,plot = F) # Suppression des silences - Semble assez chronophage
  d2        <- length(sound@left)/f # Duree finale (apr?s suppression des silences)
  
  nt    <- d_seq*f                # Nombre de points dans l'echantillon temporel
  nf    <- nt/2                   # Nombre de points dans l'echantillon frequentiel
  
  ti    <- d_seq                  # t0 pour le decoupage des morceaux (20 secondes pour les applaudissements, etc.)
  
  data  <- data.frame(matrix(nrow = (d2-ti)%/%overlap, ncol=1+1+1+1+1+3*length(labels)+1))  # Source, plage, d1, d2, meanS, Q05S, Q95S, Instrument
  
  j     <- 1
  if(d2<1e9) { # Pour pouvoir eventuellement travailler uniquement sur les "petits" morceaux
    while ((ti+d_seq)<d2) { 
      print(paste(ti, ti + d_seq, d2, dim(data_final)[1]+j-1))    
      
      sound_c     <- cutw(sound, from = ti, to = ti + d_seq, f = f)                    # Extraction d'une sequence de duree d_seq
      sound_c_f   <- as.data.frame(seewave::spec(sound_c, f = f, plot = F, fftw = T))  # D?composition en frequences et stockage dans un df
      sound_c_f$Registre <- unlist(apply(data.frame(sound_c_f$x*1e3), FUN=Registres, MARGIN=1, freqs, labels))
      
      if(max(abs(sound_c))>0) {
        out <- sound_c_f %>% 
          group_by(Registre) %>% 
          summarise(mean=mean(y), q05=quantile(y, 0.05), q95=quantile(y, 0.95))
        
        data[j, 1]                                              <- morceau
        data[j, 2]                                              <- ti
        data[j, 3:4]                                            <- c(d1, d2)
        data[j, 5:(5+length(labels)-1)]                         <- out[1,-1]
        data[j, (5+length(labels)):(5+2*length(labels)-1)]      <- out[2,-1]
        data[j, (5+2*length(labels)):(5+3*length(labels)-1)]    <- out[3,-1]
        data[j, 5+3*length(labels)]                             <- InstrInput
        data[j, 5+3*length(labels)+1]                           <- ti  
        
        j <- j+1
        
      }
      
      ti <- ti + overlap
      
      data <- as.data.frame(data)
      names(data) <- c("Morceau",
                       "t0",
                       "Duree_initiale",
                       "Duree_finale",
                       paste("A_moy", labels, sep="_"),
                       paste("A_q05", labels, sep="_"),
                       paste("A_q95", labels, sep="_"),
                       "Instrument",
                       "Temps")
    }
  }
  return(na.omit(rbind(data_final, data)))
}

NameMod <- function(char) { str_sub(char, start = 1, end = str_length(char)-1) }

NumMod  <- function(char) { str_sub(char, start = str_length(char), end = str_length(char)) }

plot_spectro <- function(sound, save, width, height, disp=F, from = from, to = to, freqs = freqs, labels = labels) {
  require(cowplot)
  
  sound_c <- cutw(sound, from = from, to = to)
  
  v <- ggspectro(sound_c, ovlp = 10, f = sound@samp.rate) + 
    geom_tile(aes(fill = amplitude)) + 
    stat_contour(geom="polygon", aes(fill=..level..), bins=30) + 
    geom_segment(aes(x = 0, y = freqs[1]/1e3, xend = to-from, yend = freqs[1]/1e3), size = 1) +
    geom_segment(aes(x = 0, y = freqs[2]/1e3, xend = to-from, yend = freqs[2]/1e3), size = 1) +
    geom_segment(aes(x = 0, y = freqs[3]/1e3, xend = to-from, yend = freqs[3]/1e3), size = 1) +
    geom_segment(aes(x = 0, y = freqs[4]/1e3, xend = to-from, yend = freqs[4]/1e3), size = 1) +
    geom_segment(aes(x = 0, y = freqs[5]/1e3, xend = to-from, yend = freqs[5]/1e3), size = 1) +
    geom_segment(aes(x = 0, y = freqs[6]/1e3, xend = to-from, yend = freqs[6]/1e3), size = 1) +
    annotate("text", x = 1, y = mean(c(10      , freqs[1]))/1e3, label = labels[1], size = 7) +
    annotate("text", x = 1, y = mean(c(freqs[1], freqs[2]))/1e3, label = labels[2], size = 7) +
    annotate("text", x = 1, y = mean(c(freqs[2], freqs[3]))/1e3, label = labels[3], size = 7) +
    annotate("text", x = 1, y = mean(c(freqs[3], freqs[4]))/1e3, label = labels[4], size = 7) +
    annotate("text", x = 1, y = mean(c(freqs[4], freqs[5]))/1e3, label = labels[5], size = 7) +
    annotate("text", x = 1, y = mean(c(freqs[5], freqs[6]))/1e3, label = labels[6], size = 7) +
    scale_fill_gradientn(name="Amplitude (dB)\n", limits=c(-30,0),
                         na.value="transparent", colours = spectro.colors(30)) +
    scale_y_log10(limits = c(NA,10)) +
    theme_bw() +
    theme(axis.text.x     = element_blank(),
          axis.text.y     = element_text(size = 16),
          axis.title.x    = element_blank(),
          axis.title.y    = element_text(size = 16),
          legend.position = "top",
          legend.text     = element_text(size = 16),
          legend.title    = element_text(size = 16))
  
  tmp <- data.frame(t=1:length(sound@left),
                    s=sound@left)
  tmp <- tmp[tmp$t%%100==0,]
  t <- ggplot(tmp) +
    geom_line(aes(x=t/sound@samp.rate, y=s/max(s))) +
    xlim(c(from, to)) +
    xlab("Time (s)") +
    ylab("Amplitude") +
    theme_shiny
  t
  final <- plot_grid(v, t, align = "hv", axis = "lr", rel_heights = c(1, 0.4), nrow = 2)
  if(disp) final
  ggsave(final, filename = save, width = width, height = height, units = "cm")
}

theme_shiny <- theme_bw() +
  theme(axis.text       = element_text(size = 16),
        axis.title.x    = element_text(size = 16),
        axis.title.y    = element_text(size = 16),
        legend.position = "top",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 16))

stack_model <- function(data, mod_stack = "mod_stack2") {
  
  mod_stack1 <- readRDS("./Modeles/Stack1.RDS")
  mod_stack2 <- readRDS("./Modeles/Stack2.RDS")
  
  don.X <- model.matrix(Instrument~., data=data)[,-1]
  don.Y <- data$Instrument
  
  mod1 <- readRDS("./Modeles/BestModel_CV.RDS")
  mod2 <- readRDS("./Modeles/GradientBoosting.RDS")
  mod3 <- readRDS("./Modeles/SVMRad_Stack.RDS")
  mod4 <- readRDS("./Modeles/Perceptron.RDS")
  
  pred1 <- as.data.frame(predict(mod1, don.X, type="prob"))
  best.iter <- gbm.perf(mod2,method="cv", plot.it=FALSE)
  pred2 <- as.data.frame(predict(mod2,data, n.trees=best.iter, type = "response"))
  pred3 <- as.data.frame(attr(predict(mod3, data, probability = T), "probabilities"))
  pred4 <- as.data.frame(predict(mod4, don.X, type = "raw"))
  
  names(pred1) <- paste0(names(pred1), "RandomForest")
  names(pred2) <- paste0(names(pred2), "GradientBoosting")
  names(pred3) <- paste0(names(pred3), "SVMradial")
  names(pred4) <- paste0(names(pred4), "Perceptron")
  
  don_stack   <- cbind(pred1, pred2, pred3, pred4, Instrument = don.Y)
  don_stack.X <- model.matrix(Instrument~., data=don_stack)[,-1]
  don_stack.Y <- don_stack$Instrument
  prev1 <- as.vector(predict(mod_stack1, newx=don_stack.X, type="class"))
  prev2 <- as.vector(predict(mod_stack2, don_stack, type="class"))
  
  if (mod_stack == "mod_stack1") {
    return(prev1)
  } else {
    return(prev2)
  }
  
}

perf_test <- function() {
  don_test   <- readRDS("./Data/don_test.RDS")
  don_test.X <- model.matrix(Instrument~., data=don_test)[,-1]
  don_test.Y <- don_test$Instrument
  
  PREV_test <- data.frame("Instrument"         = don_test.Y, 
                          "Ridge0"             = NA, 
                          "LASSO0"             = NA, 
                          "ElasticNet0"        = NA, 
                          "RandomForest1"      = NA, 
                          "RandomForest2"      = NA, 
                          "RandomForest3"      = NA, 
                          "RandomForest4"      = NA, 
                          "RandomForest5"      = NA, 
                          "RandomForest6"      = NA, 
                          "RandomForest7"      = NA, 
                          "RandomForest8"      = NA, 
                          "RandomForest9"      = NA,
                          "SVMlinear1"         = NA,
                          "SVMlinear2"         = NA,
                          "SVMlinear3"         = NA,
                          "SVMpolynomial1"     = NA,
                          "SVMpolynomial2"     = NA,
                          "SVMpolynomial3"     = NA, 
                          "SVMradial1"         = NA, 
                          "SVMradial2"         = NA, 
                          "SVMradial3"         = NA, 
                          "SVMsigmoid1"        = NA, 
                          "SVMsigmoid2"        = NA, 
                          "SVMsigmoid3"        = NA,
                          "GradientBoosting0"  = NA,
                          "Perceptron0"        = NA,
                          "StackLASSO"         = NA,
                          "StackRF"            = NA)
  
  ###################################
  ridge <- readRDS("./Modeles/Ridge0.RDS")
  PREV_test[, "Ridge0"] <- as.vector(predict(ridge, newx=don_test.X, type="class")) # Par defaut on utilise lambda.1se
  print("   RIDGE IS DONE!!")
  
  ###################################
  lasso <- readRDS("./Modeles/LASSO0.RDS")
  PREV_test[, "LASSO0"] <- as.vector(predict(lasso, newx=don_test.X, type="class")) # Par defaut on utilise lambda.1se
  print("   LASSO IS DONE!!")
  
  ###################################
  elast <- readRDS("./Modeles/ElasticNet0.RDS")
  PREV_test[, "ElasticNet0"] <- as.vector(predict(elast, newx=don_test.X, type="class")) # Par defaut on utilise lambda.1se
  print("   ELASTIC-NET IS DONE!!")
  
  for (j in 1:dim(p_rf)[1]) {
    ###################################
    name <- paste0("RandomForest",j)
    random.forest <- readRDS(paste0("./Modeles/",name,".RDS"))
    PREV_test[, name] <- as.vector(predict(random.forest, don_test.X, type="class"))
    print(paste0("   RANDOM FOREST (NTREE = ", p_rf$ntree[j], ", MTRY = ",p_rf$mtry[j], ") IS DONE!!"))
  }
  
  for (j in 1:3) {
    ###################################
    name <- paste0("SVMlinear",j)
    svm_linear <- readRDS(paste0("./Modeles/",name,".RDS"))
    PREV_test[, name] <- as.vector(predict(svm_linear, don_test, type="class"))
    print(paste0("   SVM LINEAR (C = ", cost[j], ") IS DONE!!"))
    
    ###################################
    name <- paste0("SVMpolynomial",j)
    svm_polynomial <- readRDS(paste0("./Modeles/",name,".RDS"))
    PREV_test[, name] <- as.vector(predict(svm_polynomial, don_test, type="class"))
    print(paste0("   SVM POLYNOMIAL (C = ", cost[j], ") IS DONE!!"))
    
    ###################################
    name <- paste0("SVMradial",j)
    svm_radial <- readRDS(paste0("./Modeles/",name,".RDS"))
    PREV_test[, name] <- as.vector(predict(svm_radial, don_test, type="class"))
    print(paste0("   SVM RADIAL (C = ", cost[j], ") IS DONE!!"))
    
    ###################################
    name <- paste0("SVMsigmoid",j)
    svm_sigmoid <- readRDS(paste0("./Modeles/",name,".RDS"))
    PREV_test[, name] <- as.vector(predict(svm_sigmoid, don_test, type="class"))
    print(paste0("   SVM SIGMOID (C = ", cost[j], ") IS DONE!!"))
  }
  
  ###################################
  boost <- readRDS("./Modeles/GradientBoosting.RDS")
  best.iter <- gbm.perf(boost,method="cv", plot.it=FALSE)
  tmp <- predict(boost, don_test, n.trees=best.iter, type = "response")
  PREV_test[,"GradientBoosting0"] <- levels(don_test$Instrument)[apply(tmp, FUN=which.max, MARGIN=1)]
  print("   GRADIENT BOOSTING IS DONE!!")
  
  ###################################
  if (!interactions & !carres) {
    resN <- readRDS("./Modeles/Perceptron.RDS")
    PREV_test[,"Perceptron0"] <- predict(resN, don_test.X, type = "class")
    print("   PERCEPTRON IS DONE!!")
  } else {
    print("   TOO MANY VARIABLES FOR PERCEPTRON!!")
  }
  
  ###################################
  mod_stack1 <- readRDS("./Modeles/Stack1.RDS")
  mod_stack2 <- readRDS("./Modeles/Stack2.RDS")
  PREV_test[,"StackLASSO"] <- stack_model(data=don_test, mod_stack = "mod_stack1")
  PREV_test[,"StackRF"]    <- stack_model(data=don_test, mod_stack = "mod_stack2")
  
  print("   STACK IS DONE!!")
  
  ###################################
  PREV_test <- as.data.frame(PREV_test)
  saveRDS(PREV_test, "./Modeles/PREV_test.RDS")
  
  # Selection des meilleurs modeles sur l'ensemble de test
  summary <- apply(PREV_test[,-1], 2, PredErr, Y=PREV_test[,1])
  nbest <- which.min(summary)
  names(summary)[nbest]
  
  name <- paste0("RandomForest",5)
  random.forest <- readRDS(paste0("./Modeles/",name,".RDS"))
  
  saveRDS(random.forest, "./Modeles/BestModel_test.RDS")
  
  # Meilleur modele: on retient le modele issu de la VC
  best <- readRDS("./Modeles/BestModel_test.RDS")
  
  summary_plot <- data.frame(Model   = names(summary),
                             Failure = summary)
  summary_plot$NameMod <- NameMod(summary_plot$Model)
  summary_plot$NumMod  <- as.integer(NumMod(summary_plot$Model))
  
  summary_plot_sum <- summary_plot %>% 
    group_by(NameMod) %>% 
    summarize(Failure = min(Failure))
  
  g <- ggplot(summary_plot_sum) +
    geom_bar(aes(x    = reorder(NameMod, Failure), 
                 y    = Failure*100,
                 fill = NameMod), 
             stat="identity", 
             alpha=0.75) +
    xlab("") +
    ylab("Failure (% of the test set)") +
    scale_fill_brewer(palette="Set3") +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 16),
          axis.text.y = element_text(angle = 0 , hjust=1, size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill="white"),
          legend.position = "none")
  g
  ggsave(g, filename = "./Figures/BestModels_test.jpg", width= 20, units = "cm")
  
  g <- ggplot(summary_plot_sum) +
    geom_bar(aes(x    = reorder(NameMod, Failure), 
                 y    = 100-Failure*100,
                 fill = NameMod), 
             stat="identity", 
             alpha=0.75) +
    xlab("") +
    ylab("Success (% of the test set)") +
    scale_fill_brewer(palette="Set3") +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 16),
          axis.text.y = element_text(angle = 0 , hjust=1, size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill="white"),
          legend.position = "none")
  g
  ggsave(g, filename = "./Figures/BestModels_test_perf.jpg", width= 20, units = "cm")
  
  # Table de confusion/contingence
  cm <- t(table(PREV_test[,1],PREV_test[,nbest+1]))
  df <- sweep(cm, 2, colSums(cm), FUN = "/")
  df
  df_cm <- as.data.frame(df)
  names(df_cm) <- c("Prediction", "Reference", "Freq")
  
  g <- ggplot(df_cm) +
    geom_tile(aes(x=Reference, y=Prediction, fill=Freq), alpha=0.75) +
    geom_text(aes(x=Reference, y=Prediction, label = sprintf("%1.0f", Freq*100)), vjust = 1) +
    xlab("Reference") +
    ylab("Prediction") +
    guides(fill="none") +
    coord_equal() +
    scale_fill_gradient(low="blue", high="green") +
    theme(axis.text.x  = element_text(angle = 45, hjust=1, size = 16),
          axis.text.y  = element_text(angle = 45, hjust=1, size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill="white"))
  g
  ggsave(g, filename = "./Figures/ConfusionMatrix_TestSet_test.jpg", width= 20, units = "cm")
  
  # Alluvial representation
  df_alluvial <- df_cm
  df_alluvial$Prediction <- as.factor(df_alluvial$Prediction)
  df_alluvial$Reference  <- as.factor(df_alluvial$Reference)
  
  df <- df_alluvial
  names(df) <- c("Prediction", "Original label", "Freq")
  str(df)
  
  ll <- unique(df$`Original label`)
  grid.col <- rainbow(length(ll))
  grid.col <- setNames(grid.col, ll)
  
  levs1 <- levels(df$`Original label`) 
  levs2 <- levels(df$Prediction)
  res1 <- unique(df$`Original label`)
  res2 <- unique(df$Prediction)
  cond1_cols <- grid.col[levs1[levs1 %in% res1]]
  cond2_cols <- grid.col[levs2[levs2 %in% res2]]
  columnCols <- c(cond1_cols, cond2_cols)
  stratCols <- c(rev(cond1_cols), rev(cond2_cols))
  
  df_expanded <- df[rep(row.names(df), round(df$Freq*dim(don_test)[1])), ]
  df_expanded <- df_expanded %>%
    mutate(id = row_number()) %>%
    pivot_longer(-c(Freq, id), names_to = "Condition", values_to = "label")
  
  # plot alluvial diagram
  threshold <- 0.0001
  q <- ggplot(df_expanded %>% filter(Freq > threshold), aes(x = Condition, stratum = label, alluvium = id, fill = label)) +
    geom_flow(width = 0) +
    scale_fill_manual(values = columnCols) +
    scale_color_manual(values = stratCols) +
    geom_stratum(width = 1/8, color = "white") +
    scale_x_discrete(expand = c(.25, .25)) +
    scale_y_continuous(breaks = NULL) +
    geom_text(aes(label = after_stat(stratum),
                  hjust = ifelse(Condition == "Original label", 1, 0),
                  x = as.numeric(factor(Condition)) + .075 * ifelse(Condition == "Original label", -1, 1),
                  color = after_stat(stratum)),
              stat = "stratum", 
              size = 8) +
    xlab("") +
    ylab(NULL) +
    theme_minimal() +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x=element_text(size=16),
      legend.position = "none")
  q
  ggsave(q, filename = "./Figures/Flow_test.jpg", width = 24, units = "cm")
  
}