Registres_opt_scal <- function(fi, freqs, labels) {
  return(labels[which.max(freqs[freqs-fi*1e3<=0])])
}

Registres_opt <- function(f, freqs, labels) {
  freqs  <- c(0, freqs, 1e5)
  labels <- c(labels, "Out_of_range")
  return(unlist(sapply(f[,1], FUN = Registres_opt_scal, freqs, labels)))
}

F_affect <- function(ti, d_seq, sound, freqs, labels, InstrInput) {
  require(tidyverse, warn.conflicts = F)
  require(seewave,   warn.conflicts = F)
  require(tuneR,     warn.conflicts = F)
  
  f        <- sound@samp.rate      # Sampling duration  
  d        <- length(sound@left)/f # Sound duration
  
  # print(ti)
  sound_c            <- cutw(sound, from = ti, to = ti + d_seq, f = f)
  sound_c_f          <- as.data.frame(seewave::spec(sound_c, f = f, plot = F, fftw = T))
  sound_c_f          <- sound_c_f[sound_c_f$x<10,] # We keep only frequencies below 10 kHz
  sound_c_f$Registre <- Registres_opt(sound_c_f, freqs, labels)
  
  out <- sound_c_f %>% 
    group_by(Registre) %>% 
    summarise(mean=mean(y), q05=quantile(y, 0.05), q95=quantile(y, 0.95))
  
  data <- data.frame(ncol = 24, nrow = 1)
  
  data[1, 1]                                              <- "a"
  # data[1, 1]                                              <- morceau
  data[1, 2]                                              <- ti
  data[1, 3:4]                                            <- c(d, d)
  data[1, 5:(5+length(labels)-1)]                         <- out[1,-1]
  data[1, (5+length(labels)):(5+2*length(labels)-1)]      <- out[2,-1]
  data[1, (5+2*length(labels)):(5+3*length(labels)-1)]    <- out[3,-1]
  data[1, 5+3*length(labels)]                             <- InstrInput
  data[1, 5+3*length(labels)+1]                           <- ti  # For animation purpose
  
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
  
  return(data)
}

PrepMorceau_Cut_Overlap_Par <- function(sound, Instrument, d_seq, overlap = d_seq/2, freqs, labels, workers = 8) {
  require(parallel,   warn.conflicts = F)
  require(doParallel, warn.conflicts = F)
  require(foreach,    warn.conflicts = F)
  require(progressr,  warn.conflicts = F)
  require(future,     warn.conflicts = F)
  require(doFuture,   warn.conflicts = F)

  f        <- sound@samp.rate
  d        <- length(sound@left)/f # Duration of the input sample
  time_seq <- seq(d_seq, d-d_seq, by = overlap)
  
  registerDoFuture()
  plan(multisession, workers = workers)
  with_progress({
    p <- progressor(along = 1:length(time_seq)) # Progress bar (parallel version)
    out <- foreach (time = time_seq, .combine = rbind) %dopar% {
      p()
      F_affect(ti = time, d_seq = d_seq, sound = sound, freqs = freqs, labels = labels, InstrInput = Instrument)
    }
  })
  
  # cl <- makeCluster(detectCores())
  # out <- foreach (time = time_seq, .combine = rbind) %dopar% {
  #   F_affect(ti = time, d_seq = d_seq, sound = sound, freqs = freqs, labels = labels, InstrInput = Instrument)
  # }
  # stopCluster(cl)
  
  return(out)
}

bloc_VC <- function(don, bloc, i, p_rf, cost, PREV) {
  bloc_i <- unique(bloc)[i]
  PREV_tmp <- PREV[bloc == bloc_i, ]
  print(paste0(bloc_i, " (", i, "/", nbloc, ")"))
  donA <- don[bloc!=bloc_i,]
  donT <- don[bloc==bloc_i,]
  donA.X <- model.matrix(Instrument~., data=donA)[,-1]
  donA.Y <- donA$Instrument
  donT.X <- model.matrix(Instrument~., data=donT)[,-1]
  donT.Y <- donT$Instrument
  
  ##################################
  ridge <- cv.glmnet(donA.X, donA.Y, alpha = 0, family="multinomial")
  PREV_tmp[, "Ridge0"] <- predict(ridge, newx=donT.X, type="class") # Par defaut on utilise lambda.1se
  lambdaR[i] <- ridge$lambda.1se
  print("   RIDGE IS DONE!!")
  
  ###################################
  lasso <- cv.glmnet(donA.X, donA.Y, alpha = 1, family="multinomial")
  PREV_tmp[, "LASSO0"] <- predict(lasso, newx=donT.X, type="class") # Par defaut on utilise lambda.1se
  lambdaL[i] <- lasso$lambda.1se
  print("   LASSO IS DONE!!")
  
  ###################################
  elast <- cv.glmnet(donA.X, donA.Y, alpha = 0.5, family="multinomial")
  PREV_tmp[, "ElasticNet0"] <- predict(elast, newx=donT.X, type="class") # Par defaut on utilise lambda.1se
  lambdaE[i] <- elast$lambda.1se
  print("   ELASTIC-NET IS DONE!!")
  
  for (j in 1:dim(p_rf)[1]) {
    name <- paste0("RandomForest",j)
    ###################################
    random.forest <- randomForest(Instrument~., 
                                  data  = donA,
                                  ntree = p_rf$ntree[j],
                                  mtry  = p_rf$mtry[j])
    PREV_tmp[, name] <- as.vector(predict(random.forest, donT, type="class"))
    print(paste0("   RANDOM FOREST (NTREE = ", p_rf$ntree[j], ", MTRY = ",p_rf$mtry[j], ") IS DONE!!"))
    
  }
  
  for (j in 1:3) {
    name <- paste0("SVMlinear",j)
    ###################################
    svm_linear <- svm(Instrument~.,
                      data   = donA,
                      kernel = "linear",
                      cost   = cost[j])
    PREV_tmp[, name] <- as.vector(predict(svm_linear, donT, type="class"))
    print(paste0("   SVM LINEAR (C = ", cost[j], ") IS DONE!!"))
    
    name <- paste0("SVMpolynomial",j)
    ###################################
    svm_polynomial <- svm(Instrument~.,
                          data   = donA,
                          kernel = "polynomial",
                          cost   = cost[j])
    PREV_tmp[, name] <- as.vector(predict(svm_polynomial, donT, type="class"))
    print(paste0("   SVM POLYNOMIAL (C = ", cost[j], ") IS DONE!!"))
    
    name <- paste0("SVMradial",j)
    ###################################
    svm_radial <- svm(Instrument~.,
                      data   = donA,
                      kernel = "radial",
                      cost   = cost[j])
    PREV_tmp[, name] <- as.vector(predict(svm_radial, donT, type="class"))
    print(paste0("   SVM RADIAL (C = ", cost[j], ") IS DONE!!"))
    
    name <- paste0("SVMsigmoid",j)
    ###################################
    svm_sigmoid <- svm(Instrument~.,
                       data   = donA,
                       kernel = "sigmoid",
                       cost   = cost[j])
    PREV_tmp[, name] <- as.vector(predict(svm_sigmoid, donT, type="class"))
    print(paste0("   SVM SIGMOID (C = ", cost[j], ") IS DONE!!"))
  }
  
  ###################################
  boost <- gbm(Instrument~.,
               data              = donA,
               distribution      = "multinomial",
               n.trees           = 1000,
               shrinkage         = 0.1,
               interaction.depth = 3,
               n.minobsinnode    = 10,
               cv.folds          = 5,
               keep.data         = TRUE)
  best.iter <- gbm.perf(boost,method="cv", plot.it=FALSE)
  tmp <- predict(boost, donT, n.trees=best.iter, type = "response")
  PREV_tmp[,"GradientBoosting0"] <- levels(don$Instrument)[apply(tmp, FUN=which.max, MARGIN=1)]
  print("   GRADIENT BOOSTING IS DONE!!")
  
  ###################################
  if (!interactions & !carres) {
    resN <- nnet(Instrument~.,
                 data   = donA,
                 size   = 10,
                 maxit  = 500,
                 linout = TRUE)
    PREV_tmp[,"Perceptron0"] <- predict(resN, donT, type = "class")
    print("   PERCEPTRON IS DONE!!")
  } else {
    #print("   TOO MANY VARIABLES FOR PERCEPTRON!!")
  }
  return(PREV_tmp)
}
