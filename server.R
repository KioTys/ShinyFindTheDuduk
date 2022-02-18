


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$urlLink <- renderUI({
    if(input$RegOrCtry == 'libre'){
      textInput(width = '100%' , "lienVideo", label = NULL, placeholder = "Paste your URL here")
      }else{
        selectInput("lienVideo", label = "Choose the video",
                    choices = c("Anthemes - Pierre Boulez - Diego Tosi, violon solo" = "v=W1uyLWjP5Q8",
                                "Comfortably Numb Solo - Pink Floyd - Acoustic Guitar Cover" = "v=MRqjOuhZsoY",
                                "Duduk & Oud" = "v=N8i8NvgEQgs",
                                "Duduk solo" = "v=BmC9I7uVwwc",
                                "Jimi Hendrix On An Acoustic Guitar(Only known 2 videos RARE)" = "v=P701paKEMXs",
                                "piano jazz improvisation" = "v=QBzHqW4V3lA",
                                "Simon Frick - Violin plays 'Smells like Teen Spirit' by Nirvana" = "v=h3k4ou_9cUE",
                                "Solo Duduk - Hicaz Here in the Darkness" = "v=8VpVRzSEgF4",
                                "Zigeunerweisen (Gypsy Airs) - Sarasate - Julia Fischer" = "v=FRV4UdWyESA"),
                    selected = NULL)
        }
  })

  
      
    # observeEvent(input$lienVideo, {
      observeEvent(input$submitVideo, {

        adressHtmlReactive <- eventReactive (input$submitVideo, {
        paste0('<iframe width="560" height="315" src="https://www.youtube.com/embed/',str_split(input$lienVideo, "=")[[1]][2],'?autoplay=1','"frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture", allowfullscreen></iframe>')
      })
      output$urlVideo <- renderText({
        adressHtmlReactive()
      })
      
      if(input$RegOrCtry == 'libre'){
        download_mp_from_youtube(input$lienVideo, best)
      }
      
      c <- readRDS(paste0("../Data/", str_split(input$lienVideo, "=")[[1]][2], ".RDS"))
      
    prev_pct_for_valueBox <<- c[[2]]
    morceau <<- c[[1]]


    observeEvent(input$submitVideo,{
      output$Duduk <- renderValueBox({

        valueBox(
          value = formatC(percent(prev_pct_for_valueBox$Duduk), digit = 3, format = "f"),
          subtitle = "% Duduk",
          icon = icon("check"),
          # color = if (downloadRate >= input$rateThreshold) "yellow" else "aqua"
          color = if (prev_pct_for_valueBox$Duduk > 0.7) "green" else "orange"
        )
      })
      
      output$Guitare <- renderValueBox({
        valueBox(
          value = formatC(percent(prev_pct_for_valueBox$Guitare), digit = 3, format = "f"),
          subtitle = "% Guitare",
          icon = icon("guitar"),
          color = if (prev_pct_for_valueBox$Guitare > 0.7) "green" else "orange"
        )
      })
      output$Piano <- renderValueBox({
        valueBox(
          value = formatC(percent(prev_pct_for_valueBox$Piano), digit = 3, format = "f"),
          subtitle = "% Piano",
          icon = icon("ethernet"),
          color = if (prev_pct_for_valueBox$Piano > 0.7) "green" else "orange"
        )
      })
      output$Violon <- renderValueBox({
        valueBox(
          value = formatC(percent(prev_pct_for_valueBox$Violon), digit = 3, format = "f"),
          subtitle = "% Violon",
          icon = icon("music"),
          color = if (prev_pct_for_valueBox$Violon > 0.7) "green" else "orange"
        )
      })
      

    output$plotSpectrogramme <- renderImage({
      list(src = paste0("../Figures/Spectro_", str_split(input$lienVideo, "=")[[1]][2],".jpg"),
                width = 600,
           height = 400)
    },deleteFile=FALSE
    )
    
    # L'animation : DONE !
    output$legif <- renderImage({
      list(src = paste0("../Figures/Animation_", str_split(input$lienVideo, "=")[[1]][2],".gif"),
                contentType = 'image/gif',
           width = 600,
           height = 400)
    },deleteFile=FALSE
    )
    
    ################################
    #           PAGE 3             #
    ################################
    
    # PAGE 3 : Graphes 2 par 2 (comparaison des modèles sur un morceau donné)
    output$graphe_comp <- renderPlot({graphe_predictionpar_modele(morceau)})
    
    
    }
    )
  })    
    ################################
    #           PAGE 2             #
    ################################

      # Préparation des données "summary_plot" pour le 1er graphique
      summary <- apply(PREV_test[,-1], 2, PredErr, Y=PREV_test[,1])
      nbest <- which.min(summary)
      names(summary)[nbest]
      
      summary_plot <- data.frame(Model   = names(summary),
                                 Failure = summary)
      summary_plot$NameMod <- NameMod(summary_plot$Model)
      summary_plot$NumMod  <- as.integer(NumMod(summary_plot$Model))
      
      output$plot1 <- renderPlot({graph_model_change_color(summary_plot, input$nbest_to_save)})

      output$mc_vc <- renderPlot({graph_model_matrice_confusion(PREV_test, model_selected=input$chooseModel)})
      output$mc_test <- renderPlot({graph_model_matrice_confusion(PREV_test, model_selected=input$chooseModel)})

      output$alluvial_representation <- renderPlot({graph_model_alluvial_representation(PREV_test, model_selected=input$chooseModel, don)})
      
      # Calcul des bien classés
      
      output$pct_bien_classes_vc <- renderValueBox({
        bien_classes <- calcul_bien_classes(PREV, input$chooseModel)
        valueBox(
          value = formatC(percent(bien_classes), digit = 3, format = "f"),
          subtitle = "Données Validation : % bien classés",
          icon = icon("meh-blank"),
          color = if (bien_classes > 0.8) "green" else "orange"
        )
      })
      output$pct_bien_classes_test <- renderValueBox({
        bien_classes_test <- calcul_bien_classes(PREV_test, input$chooseModel)
        valueBox(
          value = formatC(percent(bien_classes_test), digit = 3, format = "f"),
          subtitle = "Données TEST : % bien classés",
          icon = icon("thumbs-up"),
          color = if (bien_classes_test > 0.8) "green" else "orange"
        )
      })
}
