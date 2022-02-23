library(tidyverse)
# 
# theme_shiny <- theme_bw() +
#   theme(axis.text       = element_text(size = 16),
#         axis.title.x    = element_text(size = 16),
#         axis.title.y    = element_text(size = 16),
#         legend.position = "top",
#         legend.title    = element_blank(),
#         legend.text     = element_text(size = 16))


################################
#           PAGE 1             #
################################

download_mp_from_youtube <- function(link_html, best) {
  print("link_html")
  print(link_html)
  
  withProgress(message = 'MP3 se charge',
               detail = 'Cela peut prendre un moment...', value = 0, {
                 incProgress(1/5)
                 Sys.sleep(0.25)
                 system(paste0("youtube-dl ", link_html, " -x --audio-format mp3"))
               })
  withProgress(message = 'MP3 ok !',
               detail = 'Sauvegarde en cours...', value = 1/5, {
                 incProgress(1/5)
                 Sys.sleep(0.25)
                 
                 lien <- str_split(link_html, "=")[[1]][2]
                 liste <- list.files("./")
                 print(liste)
                 morceau <- paste0("./Tests/", lien, ".mp3")
                 
                 NouveauMorceau <- liste[grep(liste, pattern=".*mp3.*")]
                 
                 file.rename(from = NouveauMorceau,
                             to = paste0(lien, ".mp3"))
                 
                 file.copy(paste0("./", lien, ".mp3"), to = paste0("./Tests/", lien, ".mp3"))
                 # file.remove(paste0("./", paste0(lien, ".mp3")))
                 file.remove(paste0("./", lien, ".mp3"))
               })
  
  withProgress(message = 'Save OK',
               detail = 'Preparation des variables explicatives', value = 2/5, {
                 incProgress(1/5)
                 Sys.sleep(0.25)
                 b <- prepare_morceau_si_new(morceau, best)
                 saveRDS(b, paste0("./Data/", lien, ".RDS"))
               })
  
  withProgress(message = 'Variables OK',
               detail = 'Graphe gif en cours...', value = 3/5, {
                 incProgress(1/5)
                 Sys.sleep(0.25)
                 
                 print("Preparation de l animation")
                 graphe_dynamique_predictions(best, lien, b[[1]], chemin = "./Figures/Animation_", morceau)
               })
  withProgress(message = 'Graphe gif OK',
               detail = 'Preparation du spectrogramme', value = 4/5, {
                 incProgress(1/5)
                 Sys.sleep(0.25)
                 
                 print("Preparation du spectrogramme")
                 graphe_spectrogramme(morceau, lien, chemin = "./Figures/Spectro_")
               })
  
}



#  PAGE 1 : Graphe des fréquences
graphe_spectrogramme <- function(morceau, lien, chemin) {
  sound   <- readMP3(morceau)
  t1 <- round(length(sound@left)/sound@samp.rate/2,0)
  t2 <- min(t1+10, round(length(sound@left)/sound@samp.rate,0))
  sound_c <- cutw(sound, from = t1, to = t2, f = sound@samp.rate)
  jpeg(paste0(chemin, lien,".jpg"), width = 20, height = 15, units = "cm", res = 144)
  spectro(sound_c, sound@samp.rate, 
          flim = c(0,5), 
          flog = T, 
          fastdisp = T, 
          osc = T, 
          grid = T,
          tlab = "Temps (s)",
          flab = "Frequence (kHz)",
          title = "Spectrogramme du morceau")
  dev.off()
}


# prepare_morceau_si_new <- function(lien, best) {
prepare_morceau_si_new <- function(morceau, best) {
  instrument <- "Unknown"
  sound <- readMP3(morceau)
  labels <- c("ExtrGrave", "Grave", "BasMed", "HautsMed", "Aigu", "ExtrAigu")
  freqs  <- c(60, 200, 640, 1625, 7000, 10000)
  donnees <- PrepMorceau_Cut_Overlap_Par(sound, instrument,
                                         d_seq = 10, 
                                         overlap = 1, 
                                         freqs = freqs, 
                                         labels = labels,
                                         workers = 2)
  donnees <- donnees[,-c(1:4)]
  
  donnees$Instrument <- as.factor(donnees$Instrument)
  
  predict <- as.data.frame(predict(best, donnees, type="class"))
  names(predict) <- "PredInstr"
  print("ok")
  # browse()
  essai <- predict %>%
    group_by (PredInstr, .drop = FALSE) %>%
    summarize(compte = n()) %>%
    mutate(part = compte/nrow(predict)) %>%
    select(-compte)
  a <- dcast(essai, . ~ PredInstr, value.var = "part")
  # %>%
  #   mutate (morceau = lien)
  
  return(list(donnees, a))
}


#  PAGE 1 : Graphe dynamique prédictions
# graphe_dynamique_predictions <- function(best, lien, donnees) {
graphe_dynamique_predictions <- function(best, lien, donnees, chemin, morceau) {
  sound <- readMP3(morceau)
  
  predict2        <- as.data.frame(predict(best, donnees, type = "prob"))
  predict2$Temps  <- donnees$Temps
  predict2_cumsum <- as.data.frame(cbind(t(apply(predict2[,-5], FUN=cumsum, MARGIN=1)),
                                         Temps = predict2$Temps))
  df              <- gather(predict2,        key = "Instrument", value = "Probabilit", -5)
  df_cumsum       <- gather(predict2_cumsum, key = "Instrument", value = "Probabilit", -5)
  
  tmp <- as.vector(df_cumsum %>% 
                     filter(Instrument != "Violon") %>% 
                     select(`Probabilit`))
  Probabilit_inf <- c(rep(0, dim(df_cumsum)[1]/4), tmp$Probabilit)
  
  df_cumsum$Proba_inf <- Probabilit_inf
  
  static <- ggplot(df_cumsum) + 
    geom_ribbon(aes(x = Temps, ymin=Proba_inf*100, ymax=Probabilit*100, fill = Instrument), alpha = 0.5) +
    theme_shiny
  static
  
  dynamic <- static + transition_reveal(Temps)
  
  
  animate(dynamic, height = 500, width = 1000, duration = length(sound@left)/sound@samp.rate, nframes = round(length(sound@left)/sound@samp.rate, 0))
  
  
  anim_save(filename = paste0(chemin, lien,".gif"))
  
}

################################
#           PAGE 2             #
################################
# PAGE 2 : Graphique donnant l'erreur par type de modèle.

graph_model_change_color <- function(summary_plot, nbest_to_save) {
  failure_threshold <- summary_plot$Failure[order(summary_plot$Failure)[nbest_to_save]]
  summary_plot <- summary_plot %>%
    mutate(Keep = Failure <= failure_threshold)
  
  gg <- ggplot(summary_plot) +
    geom_bar(aes(x    = reorder(Model, Failure), 
                 y    = Failure*100,
                 fill = Keep), 
             stat="identity", 
             alpha=0.75) +
    xlab("Model") +
    ylab("Failure (% of the test set)") +
    scale_color_manual(values = c("TRUE"  = "green",
                                  "FALSE" = "red"),
                       aesthetics = c("color", "fill")) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          panel.background = element_rect(fill="white"),
          legend.position = "none")
  return(gg)
}

# PAGE 2 : Matrice de confusion

graph_model_matrice_confusion <- function(PREV, model_selected) {
  cm <- t(table(PREV[,1], PREV[[model_selected]]))
  df <- sweep(cm, 2, colSums(cm), FUN = "/")
  df
  df_cm <- as.data.frame(df)
  names(df_cm) <- c("Prediction", "Reference", "Freq")
  
  ggplot(df_cm) +
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
}

# names(df_alluvial)

# PAGE 2 : Alluvial representation
graph_model_alluvial_representation <- function(PREV, model_selected, don_test) {
  cm <- t(table(PREV[,1], PREV[[model_selected]]))
  df <- sweep(cm, 2, colSums(cm), FUN = "/")
  df_alluvial <- as.data.frame(df)
  names(df_alluvial) <- c("Prediction", "Reference", "Freq")
  df_alluvial$Prediction <- as.factor(df_alluvial$Prediction)
  df_alluvial$Reference  <- as.factor(df_alluvial$Reference)
  df <- df_alluvial
  names(df) <- c("Prediction", "Original label", "Freq")
  
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
  threshold <- 0.01
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
  return(q)
}

# PAGE 2 : Calcul du % de bien classés
calcul_bien_classes <- function(ma_table, model_selected) {
  hello <- ma_table %>%
    select("Instrument", model_selected) %>%
    rename(a = model_selected) %>%
    mutate(a = as.factor(a),
           bien_classe = ifelse(Instrument==a, 1, 0))
  
  return(sum(hello$bien_classe)/nrow(hello))
}


# Graphique du nombre de données (avant correction) par instrument
graph_donnees <- function(don){
  den <- nrow(don)
  df <- don %>%
    group_by(Instrument) %>%
    summarise(n = n())  %>%
    arrange(desc(n)) %>%
    mutate(fraction = n / den)
  
  gg <- ggplot(df, aes(x = reorder(Instrument, -fraction), y = fraction, fill = Instrument)) +
    geom_bar(stat="identity", 
             position = "dodge")+
    scale_y_continuous(labels = scales::percent) +
    coord_flip() +
    xlab("Instrument")
  ylab("Percentage of data")
  
  return (gg)  
}

################################
#           PAGE 3             #
################################

graphe_predictionpar_modele <- function(donnees) {
  pred1 <- cbind(Temps = donnees$Temps, as.data.frame(predict(mod1, donnees[,1:18], type = "prob")), Model = "Best")
  pred2 <- cbind(Temps = donnees$Temps, as.data.frame(predict(mod2, donnees[,1:18])), Model = "Perceptron")
  pred3 <- cbind(Temps = donnees$Temps, as.data.frame(predict(mod3, donnees[,1:18], n.trees=25, type = "response")[,,1]), Model = "GradientBoosting")
  pred4 <- cbind(Temps = donnees$Temps, as.data.frame(predict(mod4, model.matrix(Instrument~., data=donnees)[,-c(1, ncol(donnees))], type = "response")), Model = "LASSO")
  names(pred4) <- names(pred3)
  
  df <- rbind(pred1, pred2, pred3, pred4) %>%
    gather(key = "Instrument", value = "Proba", 2:5)
  
  g <- ggplot(df) +
    geom_line(aes(x = Temps, y = Proba*100, color = Model), size = 2) +
    facet_wrap(.~Instrument, ncol = 2, nrow = 2, scales = "free_x") +
    xlab("Time (s)") +
    ylab("Probability (%)") +
    theme_bw()+
    theme(axis.text.x=element_text(size=16),
          axis.text.y=element_text(size=16),
          axis.title.x=element_text(size=16),
          axis.title.y=element_text(size=16),
          strip.text.x = element_text(size = 16),
          legend.position = "top",
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 16))
  # ggsave(g, filename = "CompModel.jpg", width= 20, units = "cm")
  return(g)
}
