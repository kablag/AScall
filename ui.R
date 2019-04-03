library(shiny)
library(shinyMolBio)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("AScall"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("inputFile",
                HTML("Upload <b>.rdml</b>, <b>.csv</b>, <b>.lc96p</b>, <b>.xls</b> or <b>.xlsx</b> file(s):")),
      wellPanel(
        checkboxInput("preprocessCheck", "Preprocess Curves", TRUE),
        conditionalPanel(
          "input.preprocessCheck == true",
          sliderInput("bgRange", "Background Range",
                       1, 40, c(10, 20), 1))
      ),
      wellPanel(
        selectInput("ctrlMarker", "Control Marker", choices = ""),
        # numericInput("ctrlDelta", "Control ∆", 1.5),
        numericInput("cqDelta", "Cq ∆", 2),
        numericInput("cqThr", "Cq Threshold", 30),
        numericInput("rfuThr", "RFU Threshold", 1000)),
      actionButton("recalculate", "Calculate Results"),
      width = 2
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      wellPanel(
        fluidRow(
          column(6, selectInput("showKits", "Show Kits", choices = "", multiple = TRUE)),
          column(6, selectInput("showMarkers", "Show Markers", choices = "", multiple = TRUE))
        )
      ),
      fluidRow(
        column(6, uiOutput("ampCurvesUI")),
        column(6, uiOutput("pcrPlateUI"))
      ),
      dataTableOutput("globalResultsTbl")
    )
  )
)