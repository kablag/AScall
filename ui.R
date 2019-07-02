library(shiny)
library(shinyMolBio)
library(plotly)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("AScall"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("inputFile",
                "Upload .rdml, .csv, .lc96p, .xls or .xlsx file(s):",
                accept = c(".rdml", ".csv", ".lc96p", ".xls", ".xlsx")),
      wellPanel(
        tags$div(title = "Background subtraction and Cq calculation (second derivative maximum) will be applied to the curves if checked",
                 checkboxInput("preprocessCheck", "Preprocess Curves", TRUE)
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
          column(6, selectInput("showKits", "Show Kits", choices = "", multiple = TRUE)),
          column(6, selectInput("showMarkers", "Show Markers", choices = "", multiple = TRUE))
        ),
        selectInput("showSamples", "Show Samples", choices = "", multiple = TRUE)
      ),
      tabsetPanel(
        tabPanel("Summary",
                 plotlyOutput("allelicDescrPlot"),
                 dataTableOutput("summaryTbl")),
        tabPanel("Details", fluidRow(
          column(6, uiOutput("ampCurvesUI")),
          column(6, uiOutput("pcrPlateUI"))),
          dataTableOutput("detailsTbl"))
      )
    )
  )
)