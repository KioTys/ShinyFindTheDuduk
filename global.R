library(ggplot2)
library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggalluvial)
library(dplyr)
library(plotly)
library(scales)
library(seewave)
library(tuneR)
library(ggpubr)
library(png)
library(fftw)
library(seewave)
library(randomForest)
library(gganimate)
library(gifski)
library(data.table)
library(stats)
library(glmnet)
library(randomForest)
library(leaps)
library(e1071)
library(gbm) 
library(nnet) 

source("./Fonctions.R")
source("./04 - Functions for shiny.R")
source("./Fonctions_Par.R")


theme_shiny <- theme_bw() +
  theme(axis.text       = element_text(size = 16),
        axis.title.x    = element_text(size = 16),
        axis.title.y    = element_text(size = 16),
        legend.position = "top",
        legend.title    = element_blank(),
        legend.text     = element_text(size = 16))


# Chargement des données "don"
don <- readRDS("./Data/data_final - 6 bandes.rds")
# Chargement des données "PREV" sur vc
PREV <- readRDS("./Modeles/PREV_CV.RDS")

# Chargement des données "PREV" sur le test
PREV_test <- readRDS("./Modeles/PREV_test.RDS")

# Chargement des données de test
don_test <- readRDS("./Data/don_TEST.RDS")
# Chargement du meilleur modèle
best <- readRDS("./Modeles/BestModel_CV.RDS")

mod1 <- readRDS("./Modeles/BestModel_CV.RDS")
mod2 <- readRDS("./Modeles/Perceptron.RDS")
mod3 <- readRDS("./Modeles/GradientBoosting.RDS")
mod4 <- readRDS("./Modeles/LASSO0.RDS")

print("Bonjour")
