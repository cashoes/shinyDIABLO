
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(shinyBS)
library(shinydashboard)
library(shinythemes)
library(plotly)

source('helpers.R')

# UI - header ----
header <- dashboardHeader(title = "DIABLO on TCGA")

# UI - sidebar ----
sidebar <- dashboardSidebar(
  sidebarMenu(
    # menuItem("DIABLO", tabName = "diablo"),
    # menuItem("Components", tabName = "components"),
    # menuItem("Variables", tabName = "variables"),
    menuItem(
      div(
        div(
          # edit1
          style="width:75%; display:inline-block; vertical-align: middle;",
          "BiPlot"
        ),
        div(
          # edit2
          style="display:inline-block; vertical-align: middle;",
          shinyBS::bsButton("q1", label = "", icon = icon("question"),
                            style = "info", size = "small"),
          shinyBS::bsPopover(id = "q1", title = "BiPlot",
                             content = paste0("A biplot is plot which aims to represent both the observations and variables of a matrix of multivariate data on the same plot."),
                             placement = "right",
                             trigger = "click",
                             options = list(container = "body")
          )
        )
      )
      , tabName = "biplot"),
    # menuItem("Loading Vectors", tabName = "loadingVectors"),
    # menuItem("Heatmap", tabName = "heatmap"),
    menuItem(
      div(
        div(
          # edit1
          style="width:75%; display:inline-block; vertical-align: middle;",
          "Network"
        ),
        div(
          # edit2
          style="display:inline-block; vertical-align: middle;",
          shinyBS::bsButton("q2", label = "", icon = icon("question"),
                            style = "info", size = "small"),
          shinyBS::bsPopover(id = "q2", title = "Network",
                             content = paste0("Lasso a group of nodes to perform geneset enrichment analysis"),
                             placement = "right",
                             trigger = "click",
                             options = list(container = "body")
          )
        )
      )
      , tabName = "network")
    # menuItem("Circos", tabName = "circos")
  )
)

# UI - body ----
body <- dashboardBody(
  tabItems(
    tabItem("diablo",
            fixedRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         plotOutput("diablo1", height = 800)
                         ,

                         conditionalPanel(condition = "input.compareDiablo == true",
                                          plotOutput("diablo2", height = 800))

                     )
              ),
              column(width = 3,
                     box(width = NULL, status = "warning",
                         checkboxInput("compareDiablo", "Compare"),
                         p(class = "text-muted",
                           "Compare displays and contrasts another model."
                         )
                     )
              )
            )
    ),

    tabItem("components",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         plotlyOutput("compplot1", height= "auto")
                         ,
                         conditionalPanel(condition = "input.compareIndiv == true",
                                          plotOutput("indiv2")
                         )

                     ),
                     box(width = NULL,
                         verbatimTextOutput("hoverComp"),
                         verbatimTextOutput("clickComp"),
                         verbatimTextOutput("brushComp"))
              ),
              column(width = 3,
                     box(width = NULL, status = "warning",
                         checkboxInput("compareIndiv", "Compare"),
                         p(class = "text-muted",
                           "Compare displays and contrasts another model."
                         )),
                     box(width = NULL, status = "warning",
                         selectInput("selectDataComp", label = h3("Select data"),
                                     choices = dataNames,
                                     selected = 1),
                         br(),
                         radioButtons("compXComp", label = h3("X Component"),
                                      choices = as.list(1:nComp),
                                      selected = 1),
                         radioButtons("compYComp", label = h3("Y Component"),
                                      choices = as.list(1:nComp),
                                      selected = 2),
                         br(),
                         checkboxInput("showIndNames", label = "Show Ind. Names", value = FALSE),
                         p(class = "text-muted",
                           "Show individual names"
                         )
                     )
              )
            )
    ),

    tabItem("variables",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         plotlyOutput("var1", height = 800)
                         ,

                         conditionalPanel(condition = "input.compareVar == true",
                                          plotOutput("var2", height = 800))

                     ),
                     box(width = NULL,
                         DT::dataTableOutput("varTable")
                     )
              ),
              column(width = 3,
                     box(width = NULL, status = "warning",
                         checkboxInput("compareVar", "Compare"),
                         p(class = "text-muted",
                           "Compare displays and contrasts another model."
                         )),
                     box(width = NULL, status = "warning",
                         selectInput("selectDataVar", label = h3("Select data"),
                                     choices = dataNames,
                                     selected = 1),
                         br(),
                         selectInput("selectComp", label = h3("Select component"),
                                     choices = as.list(1:nComp),
                                     selected = 1),
                         p(
                           class = "text-muted",
                           paste("Select data to display in table."
                           )
                         ),
                         br(),
                         checkboxInput("showVarNames", label = "Show Var. Names", value = FALSE),
                         p(class = "text-muted",
                           "Show variable names"
                         )
                     )
              )
            )
    ),

    tabItem("biplot",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         withSpinner(plotlyOutput("biplot1", height = 800)),
                         conditionalPanel(condition = "input.compareIndiv == true",
                                          plotOutput("biplot2", height = 800))

                     )
              ),
              column(width = 3,
                     # box(width = NULL, status = "warning",
                     #     checkboxInput("compareIndiv", "Compare"),
                     #     p(class = "text-muted",
                     #       "Compare displays and contrasts another model."
                     #     )),
                     box(width = NULL, status = "warning",
                         selectInput("selectDataBi", label = h3("Select data"),
                                     choices = dataNames,
                                     selected = 1),
                         br(),
                         radioButtons("compXBi", label = h3("X Component"),
                                      choices = as.list(1:nComp),
                                      selected = 1),
                         radioButtons("compYBi", label = h3("Y Component"),
                                      choices = as.list(1:nComp),
                                      selected = 2)
                         # br(),
                         # checkboxInput("showIndNames", label = "Show Ind. Names", value = FALSE),
                         # p(class = "text-muted",
                         #   "Show individual names"
                         # )
                     )
              )
            )
    ),

    tabItem("loadingVectors",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         plotOutput("loadings1", height = 800)
                         ,

                         conditionalPanel(condition = "input.compareLoading == true",
                                          plotOutput("loadings2", height = 800))

                     )
              ),
              column(width = 3,
                     box(width = NULL, status = "warning",
                         checkboxInput("compareLoading", "Compare"),
                         p(class = "text-muted",
                           "Compare displays and contrasts another model."
                         )
                     )
              )
            )
    ),

    tabItem("heatmap",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         plotOutput("heatmap1", height = 800)
                         ,

                         conditionalPanel(condition = "input.compareHeat == true",
                                          plotOutput("heatmap2", height = 800))

                     )
              ),
              column(width = 3,
                     box(width = NULL, status = "warning",
                         checkboxInput("compareHeat", "Compare"),
                         p(class = "text-muted",
                           "Compare displays and contrasts another model."
                         )
                     )
              )
            )
    ),

    tabItem("network",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         withSpinner(plotlyOutput("network", height = 800))
                         # dataTableOutput("nodes_data_from_shiny")
                     ),
                     box(width = NULL,
                         withSpinner(dataTableOutput("brushNet"))
                         # verbatimTextOutput("brushNodes")
                         # verbatimTextOutput("clickNet")
                     )
              ),
              column(width = 3,
                     # box(width = NULL, status = "warning",
                     #     checkboxInput("compare", "Compare"),
                     #     p(class = "text-muted",
                     #       "Compare displays and contrasts another model."
                     #     )
                     # ),
                     box(width = NULL, status = "warning",
                         sliderInput("threshold", label = h3("Cutoff"),
                                     min = round(quantile(corMat, 0.50), 2),
                                     max = round(quantile(corMat, 0.99), 2),
                                     value = quantile(corMat, 0.95)),
                         # plotOutput("histoSlider", height = 100),
                         p(
                           class = "text-muted",
                           paste("Note: Cutoff removes any edges with weight less than the indicated value."
                           )
                         )
                     ),
                     box(width = NULL,
                         # valueBoxOutput("hoverNetM", width = 6),
                         # valueBoxOutput("clickNetM", width = 6),
                         valueBoxOutput(width = NULL, "density"),
                         valueBoxOutput(width = NULL, "transitivity"),
                         valueBoxOutput(width = NULL, "modularity")
                     )
                     # box(width = NULL,
                     #     fluidRow(valueBoxOutput("hoverNetM", width = 12)),
                     #     fluidRow(valueBoxOutput("clickNetM", width = 12)))
              )
            )
    ),

    tabItem("circos",
            fluidRow(
              column(width = 9,
                     box(width = NULL, solidHeader = TRUE,
                         plotOutput("circos1", height = 800),

                         conditionalPanel(condition = "input.compareCircos == true",
                                          plotOutput("circos2", height = 800))
                     )
              ),
              column(width = 3,
                     box(width = NULL, status = "warning",
                         checkboxInput("compareCircos", "Compare"),
                         p(class = "text-muted",
                           "Compare displays and contrasts another model."
                         )
                     ),
                     box(width = NULL, status = "warning",
                         sliderInput("cutoff", label = h3("Cutoff"), min = 0.5,
                                     max = 1, value = 0.90),
                         p(
                           class = "text-muted",
                           paste("Note: Cutoff removes any edges with weight less than the indicated value."
                           )
                         )
                     )
              )
            )

    )
  )
)

# UI - assemble ----
ui <- dashboardPage(
  header,
  sidebar,
  body
)

shinyUI(ui)
