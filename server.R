library(RDML)
library(tidyverse)
library(plotly)

toLog <- function(logMessage) {
  cat(sprintf("%s: %s \n", date(), logMessage))
}

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  rdml <- reactive({
    req(input$inputFile)
    toLog("Loading file")
    RDML$new(input$inputFile$datapath)
  })
  
  rdmlTbl <- reactive({
    req(rdml())
    toLog("Creating description table")
    rdml()$AsTable(cq = data$cq) %>% 
      mutate(marker = str_split(target, "_", simplify = TRUE)[, 1],
             allele = str_split(target, "_", simplify = TRUE)[, 2]) %>% 
      group_by(position) %>% 
      mutate(kit = paste(marker, collapse = ", ")) %>% 
      ungroup() %>% 
      mutate(kit = kit %>% 
               as.factor() %>% 
               as.numeric() %>% 
               sprintf("kit_%02i", .))
  })
  
  observeEvent(rdmlTbl(),
               {
                 toLog("Updating Kits select")
                 updateSelectInput(session,
                                   "showKits",
                                   choices = unique(rdmlTbl()$kit),
                                   selected = unique(rdmlTbl()$kit))
               })
  
  observeEvent(input$showKits,
               {
                 req(rdmlTbl())
                 toLog("Updating Markers select")
                 markers <- rdmlTbl() %>% 
                   filter(kit %in% input$showKits) %>% 
                   .$marker %>% 
                   unique()
                 updateSelectInput(session,
                                   "showMarkers",
                                   choices = markers,
                                   selected = markers)
               })
  
  observeEvent(rdmlTbl(),
               {
                 toLog("Updating Control Marker select")
                 updateSelectInput(session,
                                   "ctrlMarker",
                                   choices = unique(rdmlTbl()$kit),
                                   selected = unique(rdmlTbl()$kit)[1])
               })
  
  fData <- reactive({
    req(rdmlTbl())
    toLog("Getting fData")
    rdml()$GetFData(rdmlTbl(), long.table = TRUE)
  })
  
  resultsTbl <- reactive({
    req(input$calcResults)
    isolate({
      rdmlTbl() %>% 
        group_by(kit, marker, allele, sample) %>% 
        mutate(meanCq = mean(cq),
               cq_f = round(cq, digits = 2),
               meanCq_f = round(meanCq, digits = 2))
    })
  })
  
  output$pcrPlateUI <- renderUI({
    req(resultsTbl())
    toLog("Creating pcrPlate")
    pcrPlateInput("pcrPlate", 
                  plateDescription = resultsTbl() %>% 
                    dplyr::rename(sampleType = sample.type) %>% 
                    mutate(mark = sprintf("<span class='%s'></span>", sampleType)),  #whisker does not support dots!
                  wellLabelTemplate = "{{sample}}",
                  onHoverWellTextTemplate = "{{position}}\n{{sample}}\nKit: {{kit}}\nMarker: {{marker}}\nAllele: {{allele}}\nCq: {{cq_f}}\nMean Cq: {{meanCq_f}}",
                  wellClassTemplate = "{{kit}} {{sampleType}}",
                  cssText = "#{{id}} td.selected-well{border: 2px solid red !important;}
                  #{{id}} .kit_01{background-color: Plum ;}
                  #{{id}} .kit_02{background-color: LightGrey ;}
                  #{{id}} .kit_03{background-color: PaleGreen ;}
                  #{{id}} .kit_04{background-color: Salmon ;}
                  #{{id}} .kit_05{background-color: DeepSkyBlue ;}
                  #{{id}} .ntc{background-image:
    radial-gradient(#000 20%, transparent 0%),
    radial-gradient(#000 20%, transparent 0%);
  background-size: 8px 8px;
  background-position: 0 0, 4px 4px;}
                  ")
  })
  
  output$ampCurvesUI <- renderUI({
    req(fData())
    toLog("Drawing curves")
    renderAmpCurves("ampCurves", 
                    ampCurves = fData(),
                    colorBy = "marker")
  })
  
  observeEvent(
    c(input$pcrPlate, input$showMarkers),
    {
      toLog("Updating curves")
      toHideCurves <-
        which(!(rdmlTbl()$position %in% input$pcrPlate) |
                !(rdmlTbl()$kit %in% input$showKits) |
                !(rdmlTbl()$marker %in% input$showMarkers)
              # |
              #   !(rdmlTbl()$target.dyeId %in% input$showDyes)
        )
      updateCurves(session,
                   "ampCurves",
                   hideCurves = toHideCurves)
    })
  
})