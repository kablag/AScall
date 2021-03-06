source("init.R")

library(shiny)
library(shinyMolBio)
library(shinyWidgets)
library(plotly)

# Define UI for application that draws a histogram
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML(".shiny-notification {
           height: 50px;
           width: 400px;
           position:fixed;
           top: calc(50% - 25px);;
           left: calc(50% - 200px);;
           }
           "
      )
    )
  ),
  # Application title
  titlePanel("AScall"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("inputFile",
                "Upload .rdml, .csv, .lc96p, .xls or .xlsx file(s):",
                accept = c(".rdml", ".csv", ".lc96p", ".xls", ".xlsx"),
                multiple = TRUE),
      actionButton("loadExmpl", "Use Sample File"),
      tags$p(),
      wellPanel(
        tags$div(title = "Background subtraction and Cq calculation (second derivative maximum) will be applied to the curves if checked",
                 checkboxInput("preprocessCheck", "Preprocess Curves", FALSE)
        ),
        conditionalPanel(
          "input.preprocessCheck == true",
          tags$div(title = "Cycle range for background subtraction (linear part of the curves before exponentional growth)",
                   sliderInput("bgRange", "Background Range",
                               1, 40, c(5, 15), 1)
          ),
          tags$div(title = "Model type for PCR curve fitting",
                   selectInput("modelType", "Model Type",
                               c("l4", "l5", "l6", "l7", "b4", "b5", "b6", "b7"),
                               "l5")
          ),
          tags$div(title = "Cq calculation method",
                   selectInput("cqMethod", "Cq Method",
                               c("cpD2", "cpD1", "maxE", "expR", "CQ", "Cy0",
                                 "Threshold"),
                               "cpD2")
          ), 
          conditionalPanel(
            "input.cqMethod == 'Threshold'",
            tags$div(title = "Threshold value for Cq calculation",
                     numericInput("thrCq", "Threshold", 100)
            )
          )
        )
      ),
      wellPanel(
        tags$div(title = "Gene to be used as a positive control",
                 selectInput("ctrlMarker", "Control Marker", choices = "")
        ),
        # numericInput("ctrlDelta", "Control ∆", 1.5),
        tags$div(title = "Maximum cycle difference for homozygosity",
                 numericInput("cqDelta", "Cq ∆", 2)
        ),
        tags$div(title = "Mark curves as negative with Cq outside region",
                 sliderInput("cqThr", "Cq Range", 1, 50,
                             c(15, 30), step = 1)
        ),
        tags$div(title = "Mark curves as negative with RFU lower threshold",
                 numericInput("rfuThr", "RFU Threshold", 300)
        )
      ),
      conditionalPanel(
        condition = "output.enableRecalculateBtn",
        tags$div(title = "Press button to calculate results",
                 actionButton("recalculate", "Calculate Results")
        )
      ),
      conditionalPanel(
        condition = "output.resultsReady",
        tags$div(title = "Press button to generate report",
                 downloadButton("genReport", "Report")
        )
      ),
      width = 2
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      wellPanel(
        fluidRow(
          column(6, pickerInput("showKits", "Show Kits", choices = "",
                                multiple = TRUE, 
                                options = list("actions-box" = TRUE))),
          column(6, pickerInput("showMarkers", "Show Markers", choices = "",
                                multiple = TRUE, 
                                options = list("actions-box" = TRUE)))
        ),
        pickerInput("showSamples", "Show Samples", choices = "",
                    multiple = TRUE, 
                    options = list("actions-box" = TRUE))
      ),
      tabsetPanel(
        id = "tabSet",
        tabPanel("User Manual",
                 includeMarkdown("README.md")),
        tabPanel("Summary",
                 # plotlyOutput("allelicDescrPlot"),
                 plotOutput("genotypesFreqPlot"),
                 DT::dataTableOutput("summaryTbl")),
        tabPanel("Details", 
                 selectInput("showFile", "File Details", choices = ""),
                 fluidRow(
                   column(6, uiOutput("ampCurvesUI")),
                   column(6, uiOutput("pcrPlateUI"))),
                 DT::dataTableOutput("detailsTbl"))
      )
    )
  )
)