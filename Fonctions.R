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

stack_model <- function(don_test, mod_stack = "mod_stack2") {
  
  mod_stack1 <- readRDS("./Modeles/Stack1.RDS")
  mod_stack2 <- readRDS("./Modeles/Stack2.RDS")
  
  don.X <- model.matrix(Instrument~., data=don_test)[,-1]
  don.Y <- don_test$Instrument
  
  mod1 <- readRDS("./Modeles/BestModel_CV.RDS")
  mod2 <- readRDS("./Modeles/GradientBoosting.RDS")
  mod3 <- svm(Instrument~., data = don_test, kernel = "radial", cost = cost[3], probability = T)
  mod4 <- readRDS("./Modeles/Perceptron.RDS")
  
  pred1 <- as.data.frame(predict(mod1, don.X, type="prob"))
  best.iter <- gbm.perf(mod2,method="cv", plot.it=FALSE)
  pred2 <- as.data.frame(predict(mod2, don_test, n.trees=best.iter, type = "response"))
  pred3 <- as.data.frame(attr(predict(mod3, don_test, probability = T), "probabilities"))
  pred4 <- as.data.frame(predict(mod4, don.X, type = "raw"))
  
  names(pred1) <- paste0(names(pred1), "RandomForest")
  names(pred2) <- paste0(names(pred2), "GradientBoosting")
  names(pred3) <- paste0(names(pred3), "SVMradial")
  names(pred4) <- paste0(names(pred4), "Perceptron")
  
  don_stack   <- cbind(pred1, pred2, pred3, pred4, Instrument = don.Y)
  don_stack.X <- model.matrix(Instrument~., data=don_stack)[,-1]
  don_stack.Y <- don_stack$Instrument
  pred1 <- predict(mod_stack1, newx=don_stack.X, type="class")
  pred2 <- as.vector(predict(mod_stack2, don_stack, type="class"))
  
  if (mod_stack == "mod_stack1") {
    return(as.vector(pred1))
  } else {
    return(pred2)
  }
  
}
