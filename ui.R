#list of packages required
list.of.packages <- c("chipPCR", "qpcR", "DT", 
                      "tidyverse", "plotly", "RColorBrewer",
                      "shiny", "RDML", "shinyMolBio",
                      "shinyWidgets")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if (length(new.packages)) {
  cat(sprintf("Installing missing packages: %s\n", paste(new.packages, collapse = ", ")))
  install.packages(new.packages, dependencies = TRUE)
}

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
      wellPanel(
        tags$div(title = "Background subtraction and Cq calculation (second derivative maximum) will be applied to the curves if checked",
                 checkboxInput("preprocessCheck", "Preprocess Curves", FALSE)
        ),
        conditionalPanel(
          "input.preprocessCheck == true",
          tags$div(title = "Cycle range for background subtraction (linear part of the curves before exponentional growth)",
                   sliderInput("bgRange", "Background Range",
                               1, 40, c(5, 15), 1)
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
        tags$div(title = "Mark curves as negative with Cq higher threshold",
                 numericInput("cqThr", "Cq Threshold", 30)
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