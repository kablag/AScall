library(RDML)
library(tidyverse)
library(plotly)
library(doParallel)

source("generics.R")

nCores <- detectCores()
cl <- makeCluster(nCores)
registerDoParallel(cl)

onStop(function() {
  toLog("Doing application cleanup")
  stopCluster(cl)
})

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  calcResults <- reactiveVal()
  
  rdml <- reactive({
    req(input$inputFile)
    toLog("Loading file")
    withProgress(message = 'Loading file', value = 0, {
      RDML$new(input$inputFile$datapath)
    })
  })
  
  preCalcTbl <- reactive({
    req(rdml())
    rdml()$AsTable() %>% 
      addMarkersKits()
  })
  
  output$enableRecalculateBtn <- reactive({
    req(preCalcTbl())
    TRUE
  })
  outputOptions(output, 'enableRecalculateBtn', suspendWhenHidden = FALSE)
  
  observeEvent(input$recalculate,
               {
                 req(rdml())
                 
                 if (input$preprocessCheck &&
                     (is.null(calcResults()) || 
                      # don't preprocess without necessity 
                      (any(calcResults()$calcParams$bgRange != input$bgRange)))
                 ) {
                   toLog("Preprocessing curves")
                   withProgress(message = 'Preprocessing Curves', value = 0, {
                     rdml()$experiment[[1]]$run[[1]]$Preprocess(input$bgRange)
                     # for (react in rdml()$experiment[[1]]$run[[1]]$react) {
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
                 
                 
                 withProgress(message = 'Calculating Results', value = 0, {
                   for (react in rdml()$experiment[[1]]$run[[1]]$react) {
                     for (dat in react$data)
                       dat$InitEndPt()
                   }
                   
                   toLog("Creating description table")
                   dTbl <- rdml()$AsTable(cq = data$cq,
                                          endPt = data$endPt) %>% 
                     addMarkersKits()
                   toLog("Getting fData")
                   
                   
                   
                   maxCq <- rdml()$experiment[[1]]$run[[1]]$react[[1]]$data[[1]]$adp$fpoints$cyc %>% 
                     tail(1)
                   
                   # Low RFU ----------------------------------------------------------------------
                   
                   dTbl <- dTbl %>%  
                     mutate(RFU_QC = ifelse(endPt < input$rfuThr, 
                                            "Low", "Ok"))
                   # Set cq = maxCycle for curves with end point RFU lower than rfuThr | is.na(cq)
                   dTbl[dTbl$RFU_QC == "Low" | is.na(dTbl$cq), "cq"] <- maxCq
                   
                   
                   # Mark AmpStatus --------------------------------------------------------
                   
                   # noAmp if: RFU_QC != OK | is higher than cqThr
                   dTbl <- dTbl %>%  
                     mutate(ampStatus_QC = ifelse(
                       RFU_QC != "Ok" | cq > input$cqThr , 
                       "NoAmp", "Ok"))
                   
                   # Replicate match check ---------------------------------------------------
                   
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
                   
                   
                   # Kit total QC ------------------------------------------------------------
                   
                   dTbl <- dTbl %>%  
                     mutate(kit_QC = 
                     {
                       if (any(noAmpNTC_QC != "Ok"))
                         "Fail"
                       else
                         "Ok"
                     })
                   
                   
                   tmpTbl <- dTbl %>%
                     filter(kit_QC == "Ok", replicateMatch_QC == "Ok") %>%
                     group_by(kit, marker, sample) %>%
                     mutate(result = paste(allele[ampStatus_QC == "Ok"] %>% unique(), collapse = ""),
                            resultZygosity =
                              sapply(result,
                                     function(x) switch(as.character(str_length(x)),
                                                        "0" = "",
                                                        "1" = "Homo", "2" = "Hetero", "Error"))
                     )
                   dTbl <- left_join(dTbl, tmpTbl)
                   
                   # combine all errors in one column
                   qcColumns <- grep("_QC", colnames(dTbl), value = TRUE)
                   dTbl$total_QC <- apply(dTbl %>%
                                            ungroup() %>%
                                            dplyr::select(ends_with("_QC")), 1,
                                          function(x) {
                                            filterFails <- x != "Ok"
                                            paste(paste(qcColumns[filterFails],
                                                        x[filterFails], sep = ": "), collapse = "; ")
                                          })
                   
                   
                   fData <- rdml()$GetFData(dTbl, long.table = TRUE)
                 })
                 
                 calcResults(list(
                   dTbl = dTbl %>% 
                     ungroup() %>% 
                     mutate(
                       cq_f = round(cq, digits = 2),
                       meanCq_f = round(meanCq, digits = 2),
                       deltaCq_f = round(deltaCq, digits = 2),
                       result = ifelse(is.na(result), "", result)
                     ),
                   fData = fData,
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
                 updateSelectInput(session,
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
                 updateSelectInput(session,
                                   "showMarkers",
                                   choices = markers,
                                   selected = markers)
               })
  
  observeEvent(preCalcTbl(),
               {
                 toLog("Updating Samples select")
                 updateSelectInput(session,
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
  
  output$pcrPlateUI <- renderUI({
    req(calcResults())
    toLog("Creating pcrPlate")
    withProgress(message = 'Creating PCR Plate', value = 0, {
      pcrPlateInput("pcrPlate", 
                    plateDescription = calcResults()$dTbl %>% 
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
                  ",
                    interactive = TRUE)
    })
  })
  
  output$ampCurvesUI <- renderUI({
    req(calcResults())
    toLog("Drawing curves")
    withProgress(message = 'Drawing Curves', value = 0, {
    renderAmpCurves("ampCurves", 
                    ampCurves = calcResults()$fData,
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
        calcResults()$dTbl %>% 
        filter(!(position %in% input$pcrPlate) |
                 !(kit %in% input$showKits) |
                 !(marker %in% input$showMarkers) |
                 !(sample %in% input$showSamples)) %>%
        .$fdata.name
        
      
        # which(!(calcResults()$dTbl$position %in% input$pcrPlate) |
        #         !(calcResults()$dTbl$kit %in% input$showKits) |
        #         !(calcResults()$dTbl$marker %in% input$showMarkers) |
        #         !(calcResults()$dTbl$sample %in% input$showSamples)
        #       # |
        #       #   !(calcResults()$dTbl$target.dyeId %in% input$showDyes)
        # )
      updateCurves(session,
                   "ampCurves",
                   hideCurves = toHideCurves)
    })
  
  output$detailsTbl <- renderDataTable({
    req(calcResults())
    toLog("Creating DetailsTbl")
    calcResults()$dTbl %>%
      filter(position %in% input$pcrPlate &
               kit %in% input$showKits &
               marker %in% input$showMarkers &
               sample %in% input$showSamples) %>% 
      ungroup() %>% 
      dplyr::select(position, marker, kit, sample, sample.type,
                    Cq = cq_f, Mean_Cq = meanCq_f, Delta_Cq = deltaCq_f,
                    result, Zygosity = resultZygosity,
                    kit_QC, total_QC)
  })
  
  output$summaryTbl <- renderDataTable({
    req(calcResults())
    toLog("Creating ResultsTbl")
    calcResults()$dTbl %>%
      filter(kit %in% input$showKits &
               marker %in% input$showMarkers &
               sample %in% input$showSamples &
               marker != input$ctrlMarker) %>% 
      ungroup() %>% 
      dplyr::select(marker, sample, result) %>% 
      dplyr::distinct() %>% 
      spread(marker, result)
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
})



# res <- dTbl %>% 
#   group_by(fdata.name) %>% 
#   mutate(
#     cq = ifelse(is.na(cq), maxCycle, cq),
#     cq_QC = ifelse(cq <= input$cqThr, 
#                    "Ok", 
#                    sprintf("Cq %.2f > %.2f", cq, input$cqThr)) ) %>% 
#   group_by(kit, marker, allele, sample) %>% 
#   mutate(
#     meanCq = mean(cq),
#     meanCq_QC = ifelse(meanCq <= input$cqThr, 
#                        "Ok", 
#                        sprintf("Mean Cq %.2f > %.2f", meanCq, input$cqThr)),
#     deltaCq = max(cq) - min(cq),
#     # check that Cq ∆ is lower than thr value (ctrl marker has own thr)
#     deltacq_QC = { 
#       ifelse(deltaCq <= ifelse(marker == input$ctrlMarker,
#                                input$cqDelta, 
#                                input$cqDelta),
#              "Ok",
#              sprintf("Mean Cq %.2f > %.2f", meanCq, ifelse(marker == input$ctrlMarker,
#                                                            input$cqDelta, 
#                                                            input$cqDelta))
#       )
#     },
#     # Preformat for output
#     cq_f = round(cq, digits = 2),
#     meanCq_f = round(meanCq, digits = 2),
#     deltaCq_f = round(deltaCq, digits = 2)) %>% 
#   group_by(position) %>% 
#   mutate(
#     ctrlMarker_QC = 
#       if (cq_QC[marker == input$ctrlMarker] != "Ok" ||
#           meanCq_QC[marker == input$ctrlMarker] != "Ok" ||
#           deltacq_QC[marker == input$ctrlMarker] != "Ok") "Bad Control Marker"
#     else "Ok"
#   ) %>%
#   group_by(kit) %>% 
#   mutate(
#     ntc_QC = if (any(cq_QC[sample.type == "ntc"] == "Ok"))
#       "Bad NTC" else "Ok",
#     kit_QC = ifelse(ntc_QC != "Ok", "Kit Error", "Ok")
#   ) %>% 
#   group_by(kit, allele) %>% 
#   mutate(
#     pos_QC = if (all(ctrlMarker_QC[sample.type == "pos"] != "Ok") &
#                  xor(all(cq_QC[sample.type == "pos" & 
#                                grepl("ref", sample) & 
#                                marker != input$ctrlMarker] == "Ok") & 
#                      all(meanCq_QC[sample.type == "pos" & 
#                                    grepl("ref", sample) & 
#                                    marker != input$ctrlMarker] == "Ok") &
#                      all(deltacq_QC[sample.type == "pos" & 
#                                     grepl("ref", sample) & 
#                                     marker != input$ctrlMarker] == "Ok"),
#                      all(cq_QC[sample.type == "pos" & 
#                                grepl("alt", sample) & 
#                                marker != input$ctrlMarker] == "Ok") & 
#                      all(meanCq_QC[sample.type == "pos" & 
#                                    grepl("alt", sample) & 
#                                    marker != input$ctrlMarker] == "Ok") &
#                      all(deltacq_QC[sample.type == "pos" & 
#                                     grepl("alt", sample) & 
#                                     marker != input$ctrlMarker] == "Ok")))
#     "Ok" else "Bad Pos"
#     # TODO: check alt = Ok & alt = Ok   or  ref = Ok & ref = Ok
#     
#   )
# fRes <- res %>% 
#   # filter(ctrlMarker_QC == "Ok", kit_QC == "Ok") %>% 
#   group_by(kit, marker, sample) %>% 
#   mutate(result = paste(allele[cq_QC == "Ok"] %>% unique(), collapse = ""),
#          resultZygosity = 
#            sapply(result, 
#                   function(res) switch(as.character(str_length(res)),
#                                        "0" = "",
#                                        "1" = "Homo", "2" = "Hetero", "Error"))
#   )
# res <- left_join(res, fRes)
# # testTbl <<- res
# # combine all errors in one column
# res$total_QC <- apply(res %>% 
#                         ungroup() %>% 
#                         dplyr::select(ends_with("_QC")), 1,
#                       function(x) paste(x[x != "Ok"], collapse = ", "))
# res