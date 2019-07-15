library(RDML)
library(tidyverse)
library(plotly)
# library(doParallel)
library(RColorBrewer)

source("generics.R")

# nCores <- detectCores()
# cl <- makeCluster(nCores)
# registerDoParallel(cl)

onStop(function() {
  toLog("Doing application cleanup")
  # stopCluster(cl)
})

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  calcResults <- reactiveVal()
  
  rdmls <- reactive({
    req(input$inputFile)
    toLog(str_pad("Loading Files",
                  40, "both", pad = c("-")))
    withProgress(message = "Loading file", value = 0, {
      rdmls <- lapply(seq_along(input$inputFile$datapath),
                      function(i) {
                        msg <- paste("Loading File:", input$inputFile$name[i])
                        toLog(msg)
                        incProgress(0,
                                    msg)
                        rdml <- RDML$new(input$inputFile$datapath[i])
                        incProgress(1/length(input$inputFile$datapath))
                        rdml
                      })
      # rdmls <- foreach(dpath = input$inputFile$datapath) %dopar% {
      #   library(RDML)
      #   source("generics.R")
      #   RDML$new(dpath)
      # }
    })
    names(rdmls) <- input$inputFile$name
    rdmls
  })
  
  observeEvent(rdmls(),
               {
                 toLog("Updating showFile")
                 updateSelectInput(session,
                                   "showFile",
                                   choices = names(rdmls()),
                                   selected = names(rdmls())[1])
               })
  
  preCalcTbl <- reactive({
    req(rdmls())
    tbl <- data.frame()
    for (fname in names(rdmls())) {
      subTbl <- rdmls()[[fname]]$AsTable() %>% 
        mutate(fname = fname) %>% 
        addMarkersKits()
      tbl <- bind_rows(tbl, subTbl)
    }
    tbl
  })
  
  output$enableRecalculateBtn <- reactive({
    req(preCalcTbl())
    TRUE
  })
  outputOptions(output, "enableRecalculateBtn", suspendWhenHidden = FALSE)
  
  observeEvent(input$recalculate,
               {
                 req(rdmls())
                 
                 if (input$preprocessCheck &&
                     (is.null(calcResults()) || 
                      # don't preprocess without necessity 
                      (any(calcResults()$calcParams$bgRange != input$bgRange)))
                 ) {
                   toLog(str_pad("Preprocessing curves",
                                 40, "both", pad = c("-")))
                   withProgress(message = 'Preprocessing Curves', value = 0, {
                     for (fname in names(rdmls())) {
                       toLog(paste("File", fname))
                       rdmls()[[fname]]$experiment[[1]]$run[[1]]$Preprocess(input$bgRange)
                       incProgress(1/length(rdmls()))
                     }
                     # for (react in rdmls()$experiment[[1]]$run[[1]]$react) {
                     #   bgRange <- input$bgRange
                     #   results <- 
                     #     # for (dat in react$data) {
                     #     # dat$Preprocess(bgRange)
                     #     foreach(dat = react$data) %dopar% {
                     #     dat$Preprocess(bgRange)
                     #   }
                     #   for (result in results[results != "Ok"]) {
                     #     toLog(result[[2]], result[[1]])
                     #   }
                     # }
                   })
                 } else {
                   toLog("Passing curves preprocess")
                 }
                 
                 globalDescriptionTbl <- data.frame()
                 globalFluoData <- list()
                 
                 withProgress(message = 'Calculating Results', value = 0, {
                   
                   for (fname in names(rdmls())) {
                     str_pad(paste("Calculating Results for:", fname),
                             80, "both", pad = c("-")) %>% 
                       toLog()
                     for (react in rdmls()[[fname]]$experiment[[1]]$run[[1]]$react) {
                       for (dat in react$data)
                         dat$InitEndPt()
                     }
                     
                     toLog("Creating description table")
                     dTbl <- rdmls()[[fname]]$AsTable(cq = data$cq,
                                                      endPt = data$endPt) %>% 
                       addMarkersKits()
                     
                     maxCq <- rdmls()[[fname]]$experiment[[1]]$run[[1]]$react[[1]]$data[[1]]$adp$fpoints$cyc %>% 
                       tail(1)
                     
                     # Low RFU ----------------------------------------------------------------------
                     toLog("Check Low RFU")
                     dTbl <- dTbl %>%  
                       mutate(RFU_QC = ifelse(endPt < input$rfuThr, 
                                              "Low", "Ok"))
                     # Set cq = maxCycle for curves with end point RFU lower than rfuThr | is.na(cq)
                     dTbl[dTbl$RFU_QC == "Low" | is.na(dTbl$cq), "cq"] <- maxCq
                     
                     
                     # Mark AmpStatus --------------------------------------------------------
                     toLog("Mark AmpStatus")
                     # noAmp if: RFU_QC != OK | is higher than cqThr
                     dTbl <- dTbl %>%  
                       mutate(ampStatus_QC = ifelse(
                         RFU_QC != "Ok" | cq > input$cqThr , 
                         "NoAmp", "Ok"))
                     
                     # Replicate match check ---------------------------------------------------
                     toLog("Replicate match check")
                     # All replicates have to be NoAmp or (Ok and Cq∆ lower than cqDelta)
                     dTbl <- dTbl %>%
                       group_by(kit, marker, allele, sample) %>% 
                       mutate(
                         meanCq = mean(cq),
                         deltaCq = max(cq) - min(cq),
                         replicateMatch_QC = { 
                           if ((all(ampStatus_QC == "Ok") && deltaCq[1] < input$cqDelta) ||
                               all(ampStatus_QC != "Ok") 
                               # !!!! temporary hack for control marker in positive sample !!!!
                               || (sample.type == "pos" && marker == input$ctrlMarker)
                           ) "Ok"
                           else "Fail"
                         })
                     
                     
                     # Kit NTC noAmp -----------------------------------------------------------
                     toLog("Kit NTC noAmp")
                     # All NTC reactions have to be noAmp 
                     dTbl <- dTbl %>%  
                       group_by(kit) %>% 
                       mutate(noAmpNTC_QC = 
                       {
                         if (any(ampStatus_QC[sample.type == "ntc"] == "Ok"))
                           "Fail"
                         else
                           "Ok"
                       })
                     
                     
                     # Control Marker QC ---------------------------------------------------------------
                     toLog("Control Marker QC")
                     dTbl <- dTbl %>%  
                       group_by(kit, sample) %>% 
                       mutate(ctrlMarker_QC = 
                       {
                         # ampOk (and replicateMatch_QC ok )for ctrlMarker 
                         if (any(replicateMatch_QC[marker == input$ctrlMarker] != "Ok") ||
                             any(ampStatus_QC[marker == input$ctrlMarker] == "NoAmp"))
                           "Fail"
                         else
                           "Ok"
                       })
                     
                     # Kit total QC ------------------------------------------------------------
                     toLog("Kit total QC")
                     dTbl <- dTbl %>%  
                       mutate(kit_QC = 
                       {
                         # noAmp for NTC 
                         if (any(noAmpNTC_QC != "Ok"))
                           "Fail"
                         else
                           "Ok"
                       })
                     
                     
                     # Results Calc ------------------------------------------------------------
                     toLog("Results Calc")
                     genResult <- function(okAlleles) {
                       okAlleles <-  unique(okAlleles)
                       if (length(okAlleles) == 1) {
                         okAlleles <- c(okAlleles, okAlleles)
                       }
                       paste(okAlleles, collapse = "")
                     }
                     
                     tmpTbl <- dTbl %>%
                       filter(kit_QC == "Ok",
                              replicateMatch_QC == "Ok",
                              ctrlMarker_QC == "Ok") %>%
                       group_by(kit, marker, sample) %>%
                       mutate(
                         result = genResult(allele[ampStatus_QC == "Ok"]),
                         resultZygosity =
                           sapply(result,
                                  function(x) switch(
                                    as.character(str_length(x)),
                                    "0" = "",
                                    "1" = "Homo",
                                    "2" = "Hetero",
                                    "Error")
                                  )
                       )
                     dTbl <- left_join(dTbl, tmpTbl)
                     dTbl$result[dTbl$result == ""] <- NA

                     # Combine All Errors in One Column ----------------------------------------
                     toLog("Combine All Errors in One Column")
                     qcColumns <- grep("_QC", colnames(dTbl), value = TRUE)
                     dTbl$total_QC <- 
                       apply(dTbl %>%
                               ungroup() %>% 
                               dplyr::select(ends_with("_QC")), 1,
                             function(x) {
                               filterFails <- x != "Ok"
                               paste(paste(qcColumns[filterFails],
                                           x[filterFails], sep = ": "),
                                     collapse = "; ")
                             })
                     
                     toLog("Getting fData")
                     fData <- rdmls()[[fname]]$GetFData(dTbl, long.table = TRUE)
                     
                     dTbl <- dTbl %>% 
                       mutate(fname = fname)
                     globalDescriptionTbl <- bind_rows(globalDescriptionTbl, dTbl)
                     globalFluoData[[fname]] <- fData
                     incProgress(1/length(rdmls()))
                   }
                 })
                 names(globalFluoData) <- names(rdmls())
                 calcResults(list(
                   dTbl = globalDescriptionTbl %>% 
                     ungroup() %>% 
                     mutate(
                       cq_f = round(cq, digits = 2),
                       meanCq_f = round(meanCq, digits = 2),
                       deltaCq_f = round(deltaCq, digits = 2),
                       result = ifelse(is.na(result), 
                                       "!NA!", 
                                       result)
                     ),
                   fData = globalFluoData,
                   calcParams = list(
                     preprocessCheck = input$preprocessCheck,
                     bgRange = input$bgRange,
                     ctrlMarker = input$ctrlMarker,
                     cqDelta = input$cqDelta,
                     cqThr = input$cqThr,
                     rfuThr = input$rfuThr
                   )
                 ))
               }
  )
  
  observeEvent(preCalcTbl(),
               {
                 toLog("Updating Kits select")
                 updatePickerInput(session,
                                   "showKits",
                                   choices = unique(preCalcTbl()$kit),
                                   selected = unique(preCalcTbl()$kit))
               })
  
  observeEvent(input$showKits,
               {
                 req(preCalcTbl())
                 toLog("Updating Markers select")
                 markers <- preCalcTbl() %>% 
                   filter(kit %in% input$showKits) %>% 
                   .$marker %>% 
                   unique()
                 updatePickerInput(session,
                                   "showMarkers",
                                   choices = markers,
                                   selected = markers)
               })
  
  observeEvent(preCalcTbl(),
               {
                 toLog("Updating Samples select")
                 updatePickerInput(session,
                                   "showSamples",
                                   choices = unique(preCalcTbl()$sample),
                                   selected = unique(preCalcTbl()$sample))
               })
  
  
  observeEvent(preCalcTbl(),
               {
                 toLog("Updating Control Marker select")
                 updateSelectInput(session,
                                   "ctrlMarker",
                                   choices = unique(preCalcTbl()$marker),
                                   selected = unique(preCalcTbl()$marker)[1])
               })
  
  
  
  # resultsTbl <- reactive({
  #   req(calcResults()) #req(input$calcResults)
  #   isolate({
  #     
  #   })
  # })
  
  
  preDetailsTbl <- reactive({
    req(calcResults())
    calcResults()$dTbl %>% 
      filter(fname == input$showFile)
  })
  
  output$pcrPlateUI <- renderUI({
    req(preDetailsTbl())
    toLog("Creating pcrPlate")
    withProgress(message = 'Creating PCR Plate', value = 0, {
      kitsNames <- unique(preDetailsTbl()$kit)
      nKits <- length(kitsNames)
      kitsColors <- brewer.pal(nKits, "Accent")[1:nKits]
      kitsColorsCSS <- paste(
        sprintf("#{{id}} .%s{background-color: %s ;}",
                kitsNames, kitsColors), collapse = "")
      legendText <- tags$div(
              lapply(kitsNames,
                     function(kname) tags$span(class = kname, kname)),
              tags$span(class = "ntc", "NTC")
              )
      
      pcrPlateInput("pcrPlate", 
                    plateDescription = preDetailsTbl() %>% 
                      filter(kit %in% input$showKits & 
                               marker %in% input$showMarkers &
                               sample %in% input$showSamples) %>% 
                      dplyr::rename(sampleType = sample.type) %>% #whisker does not support dots!
                      group_by(position) %>% 
                      mutate(#mark = sprintf("<span class='%s'></span>", sampleType),
                        onHoverTbl = 
                          data.frame(marker, 
                                     allele,
                                     cq_f,
                                     meanCq_f,
                                     deltaCq_f,
                                     total_QC,
                                     result) %>% 
                          apply(1, function(x) sprintf("Marker %s Allele %s:\nCq = %s, Mean Cq = %s, ∆ Cq = %s, Errors = {%s}",
                                                       x["marker"], x["allele"],
                                                       x["cq_f"],
                                                       x["meanCq_f"],
                                                       x["deltaCq_f"],
                                                       x["total_QC"])) %>% 
                          paste(collapse = "\n"),
                        kitError = ifelse(kit_QC == "Kit Error", "✖︎", "")),
                    wellLabelTemplate = "<b>{{kitError}}</b> {{sample}}",
                    onHoverWellTextTemplate = "{{position}}\n{{sample}}\n{{sampleType}}\n{{kit}}\n{{onHoverTbl}}",
                    wellClassTemplate = "{{kit}} {{kit_QC}} {{sampleType}}",
                    cssText = paste0("#{{id}} td.selected-well{border: 2px solid red !important;}",
                                     kitsColorsCSS,
                  "#{{id}} .ntc{background-image:
    radial-gradient(#000 20%, transparent 0%),
    radial-gradient(#000 20%, transparent 0%);
  background-size: 8px 8px;
  background-position: 0 0, 4px 4px;}"),
                  legend = legendText,
                    interactive = TRUE)
    })
  })
  
  output$ampCurvesUI <- renderUI({
    req(calcResults())
    toLog("Drawing curves")
    withProgress(message = 'Drawing Curves', value = 0, {
      renderAmpCurves("ampCurves", 
                      ampCurves = calcResults()$fData[[input$showFile]],
                      colorBy = "marker",
                      interactive = TRUE)
    })
  })
  
  observeEvent(
    c(input$pcrPlate, input$showMarkers, input$showSamples),
    {
      req(calcResults())
      toLog("Updating curves")
      toHideCurves <-
        preDetailsTbl() %>% 
        filter(!(position %in% input$pcrPlate) |
                 !(kit %in% input$showKits) |
                 !(marker %in% input$showMarkers) |
                 !(sample %in% input$showSamples)) %>%
        .$fdata.name
      
      updateCurves(session,
                   "ampCurves",
                   hideCurves = toHideCurves)
    })
  
  observeEvent(input$hoverfDataName,
               {
                 hoverfDataName <- input$hoverfDataName
                 updateCurves(session,
                              "ampCurves",
                              highlightCurves = hoverfDataName)
                 updatePcrPlateInput(
                   session,
                   "pcrPlate",
                   highlightning = str_sub(hoverfDataName, end = 3))
                 paste("Highlighted:", hoverfDataName)
               })
  
  output$detailsTbl <- DT::renderDataTable({
    req(preDetailsTbl())
    toLog("Creating DetailsTbl")
    DT::datatable(
      preDetailsTbl() %>%
        filter(position %in% input$pcrPlate &
                 kit %in% input$showKits &
                 marker %in% input$showMarkers &
                 sample %in% input$showSamples) %>% 
        ungroup() %>% 
        dplyr::select(fdata.name, position, marker, kit, sample, sample.type,
                      Cq = cq_f, Mean_Cq = meanCq_f, Delta_Cq = deltaCq_f,
                      result, Zygosity = resultZygosity,
                      kit_QC, total_QC),
      rownames = FALSE,
      options = list(
        rowCallback = DT::JS('function(row, data) {
                              $(row).mouseenter(function(){
                              Shiny.onInputChange("hoverfDataName", data[0]);
                              });
                              $(row).mouseout(function(){
                              Shiny.onInputChange("hoverfDataName", "");
                              });}'),
        #hide column fdata.name
        columnDefs = list(list(visible = FALSE, targets = c(0))) 
      )
    )
  })
  
  output$summaryTbl <- DT::renderDataTable({
    req(calcResults())
    toLog("Creating ResultsTbl")
    DT::datatable(calcResults()$dTbl %>%
                    filter(kit %in% input$showKits &
                             marker %in% input$showMarkers &
                             sample %in% input$showSamples &
                             marker != input$ctrlMarker) %>% 
                    ungroup() %>% 
                    dplyr::select(marker, sample, result) %>% 
                    dplyr::distinct() %>% 
                    spread(marker, result))
  })
  
  output$genotypesFreqPlot <- renderPlot({
    req(calcResults())
    toLog("Creating genotypesFreqPlot")
    calcResults()$dTbl %>% 
      filter(sample.type == "unkn" & #result != "" &
               kit %in% input$showKits &
               marker %in% input$showMarkers &
               sample %in% input$showSamples &
               marker != input$ctrlMarker) %>%
      group_by(marker, sample) %>%
      summarise(result = result[1]) %>% 
      ggplot(aes(x = marker)) +
      geom_bar(aes( fill = result)) +
      geom_text(aes(label = ..count.., group = result),
                stat = "count", position = position_stack(0.4),
                color = "white") +
      geom_text(aes(label = result, group = result),
                stat = "count", position = position_stack(0.6),
                color = "white") +
      theme_bw() +
      theme(legend.position = "none")
  })
  
  output$allelicDescrPlot <- renderPlotly({
    req(calcResults())
    toLog("Creating allelicDescrPlot")
    tempTbl <- calcResults()$dTbl %>%
      filter(kit %in% input$showKits &
               marker %in% input$showMarkers &
               sample %in% input$showSamples &
               marker != input$ctrlMarker) %>% 
      dplyr::select(marker, allele, sample, result, meanCq_f) %>%
      group_by(marker) %>% 
      mutate(alleleN = allele %>% as.factor() %>% as.numeric()) %>% 
      distinct() %>% 
      filter(marker != "B2m") %>% 
      unite("sample_marker_result",sample, marker, result) %>%
      dplyr::select(-allele) %>% 
      spread(key = alleleN, value = meanCq_f)
    colnames(tempTbl)[c(2,3)] <- c("Allele_1", "Allele_2")
    # rtbl3 <<- tempTbl
    # ggplot(tempTbl) +
    # geom_point(aes(x = Allele_1, y = Allele_2))
    plot_ly() %>%
      add_trace(data = tempTbl,
                x = ~Allele_1, y = ~Allele_2,
                split = ~sample_marker_result,
                text = ~sample_marker_result,
                hoverinfo = 'text',
                marker = list(color = 'blue')
      ) %>% 
      plotly::layout(showlegend = FALSE)
  })
  
  observeEvent(
    input$pcrPlate_hover,
    {
      fdataNames <- preDetailsTbl() %>%
        filter(position %in% input$pcrPlate_hover) %>%
        .$fdata.name
      updateCurves(session,
                   "ampCurves",
                   highlightCurves = fdataNames)
    }
  )
})