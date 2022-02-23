


dashboardPage(
    dashboardHeader(title = "Find The Duduk !"),
    dashboardSidebar(
        sidebarMenu(
            menuItem("Predit nouveau modele", tabName = "dashboard1", icon = icon("th")),
            menuItem("Models overview", tabName = "dashboard2", icon = icon("th")),
            menuItem("Best models in details", tabName = "dashboard3", icon = icon("th"))
        )
    ),
    dashboardBody(
        tabItems(
            tabItem(tabName = "dashboard1",
                    h2("Please paste a youtube URL, we will tell you which instruments are playing !"),

                    
                    ###
                    fluidRow(
                      radioButtons(inputId = 'RegOrCtry', 
                                   label = 'Choisir entre une video libre ou plusieurs choix', 
                                   choices = c('libre','liste')),
                      box(uiOutput("urlLink"))
                    ),
                    actionButton("submitVideo", label = "Launch video"),
                    ###
                    fluidRow(
                        box(
                        width = 6,
                        htmlOutput("urlVideo")
                        )
                    ),
                    br(),
                    fluidRow(
                      valueBoxOutput("Duduk", width = 3),
                      valueBoxOutput("Guitare", width = 3),
                      valueBoxOutput("Piano", width = 3),
                      valueBoxOutput("Violon", width = 3)),
                    br(),
                    fluidRow(
                      box(width = 6,
                             h3("Répartition des instruments par 10 secondes d'écoute"),
                             imageOutput("legif")
                             ),
                      box(width = 6,
                             h3("Graphique des fréquences"),
                             imageOutput("plotSpectrogramme")
                    )
            )),
            tabItem(tabName = "dashboard2",
                    tabsetPanel(tabPanel("Validation",
                                         h2("Les performances de nos modèles"),
                                         fluidRow(
                                           box(width = NULL,
                                               sliderInput("nbest_to_save","nb of models retenus",1,8,3),
                                               plotOutput("plot1"))
                                         )),
                                tabPanel("Test",
                                         h2("Les performances de nos modèles"),
                                         fluidRow(
                                           box(width = NULL,
                                               sliderInput("nbest_to_save","nb of models retenus",1,8,3),
                                               plotOutput("plot2"))
                                         ))),
                    h2("Résultats par modèle"),
                    fluidRow(
                      column(width = 6,
                      selectInput("chooseModel", label = "Choose the model",
                                  choices = c("Random Forest" = "RandomForest1",
                                              "SVM Radial" = "SVMradial1",
                                              "Gradient Boosting" = "GradientBoosting0",
                                              "Perceptron" = "Perceptron0",
                                              "SVM sigmoid" = "SVMsigmoid1",
                                              "SVM linear" = "SVMlinear1"))),
                      valueBoxOutput("pct_bien_classes_vc", width = 2),
                      valueBoxOutput("pct_bien_classes_test", width = 2)),
                    fluidRow(
                      column(width = 6,
                             plotOutput("mc_vc")),
                      column(width = 6,
                             plotOutput("alluvial_representation"))
                      )
                    ),
            tabItem(tabName = "dashboard3",
                    h2("Comparaison des prévisions sur le morceau choisi"),
                    fluidRow(
                      box(width = 12, height = 550, plotOutput("graphe_comp"))
                    )
            )
        )
    )
)
