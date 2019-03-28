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
      selectInput("showKits", "Show Kits", choices = "", multiple = TRUE),
      selectInput("showMarkers", "Show Markers", choices = "", multiple = TRUE),
      selectInput("ctrlMarker", "Control Marker", choices = ""),
      numericInput("ctrlDelta", "Control ∆", 1.5),
      numericInput("trgtDelta", "Target ∆", 2),
      numericInput("cqThr", "Cq Threshold", 30),
      actionButton("calcResults", "Calculate Results"),
      width = 2
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      fluidRow(
        column(6, uiOutput("ampCurvesUI")),
        column(6, uiOutput("pcrPlateUI"))
      )
    )
  )
)