library(RDML)
library(tidyverse)
library(plotly)
library(doParallel)

source("generics.R")

nCores <- detectCores()
cl <- makeCluster(nCores)
registerDoParallel(cl)

shinyServer(function(input, output, session) {
  session$onSessionEnded(function() {
    stopApp()
  })
  
  calcResults <- reactiveVal()
  
  rdml <- reactive({
    req(input$inputFile)
    toLog("Loading file")
    RDML$new(input$inputFile$datapath)
  })
  
  preCalcTbl <- reactive({
    req(rdml())
    rdml()$AsTable() %>% 
      addMarkersKits()
  })
  
  observeEvent(input$recalculate,
    {
      req(rdml())
      if (input$preprocessCheck &&
          (is.null(calcResults()) || 
           # don't preprocess without necessity 
           (calcResults()$calcParams$bgRange != input$bgRange))) {
        toLog("Preprocessing curves")
        for (react in rdml()$experiment[[1]]$run[[1]]$react) {
          bgRange <- input$bgRange
          results <- foreach(dat = react$data) %dopar% {
            dat$Preprocess(bgRange)
          }
          for (result in results[results != "Ok"]) {
            toLog(result[[2]], result[[1]])
          }
        }
      } else {
        toLog("Passing curves preprocess")
      }
      toLog("Creating description table")
      dTbl <- rdml()$AsTable(cq = data$cq) %>% 
        addMarkersKits()
      toLog("Getting fData")
      fData <- rdml()$GetFData(dTbl, long.table = TRUE)
      calcResults(list(
        dTbl = dTbl,
        fData = fData,
        calcParams = list(
          input$preprocessCheck,
          input$bgRange,
          input$ctrlMarker,
          input$cqDelta,
          input$cqThr,
          input$rfuThr
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
                 toLog("Updating Control Marker select")
                 updateSelectInput(session,
                                   "ctrlMarker",
                                   choices = unique(preCalcTbl()$marker),
                                   selected = unique(preCalcTbl()$marker)[1])
               })
  
  
  
  resultsTbl <- reactive({
    req(calcResults()) #req(input$calcResults)
    isolate({
      maxCycle <- calcResults()$fData$cyc %>% max()
      res <- calcResults()$dTbl %>% 
        group_by(fdata.name) %>% 
        mutate(
          cq = ifelse(is.na(cq), maxCycle, cq),
          cq_QC = ifelse(cq <= input$cqThr, 
                         "Ok", 
                         sprintf("Cq %.2f > %.2f", cq, input$cqThr)) ) %>% 
        group_by(kit, marker, allele, sample) %>% 
        mutate(
          meanCq = mean(cq),
          meanCq_QC = ifelse(meanCq <= input$cqThr, 
                             "Ok", 
                             sprintf("Mean Cq %.2f > %.2f", meanCq, input$cqThr)),
          deltaCq = max(cq) - min(cq),
          # check that Cq ∆ is lower than thr value (ctrl marker has own thr)
          deltacq_QC = { 
            ifelse(deltaCq <= ifelse(marker == input$ctrlMarker,
                                     input$cqDelta, 
                                     input$cqDelta),
                   "Ok",
                   sprintf("Mean Cq %.2f > %.2f", meanCq, ifelse(marker == input$ctrlMarker,
                                                                 input$cqDelta, 
                                                                 input$cqDelta))
            )
          },
          # Preformat for output
          cq_f = round(cq, digits = 2),
          meanCq_f = round(meanCq, digits = 2),
          deltaCq_f = round(deltaCq, digits = 2)) %>% 
        group_by(position) %>% 
        mutate(
          ctrlMarker_QC = 
            if (cq_QC[marker == input$ctrlMarker] != "Ok" ||
                meanCq_QC[marker == input$ctrlMarker] != "Ok" ||
                deltacq_QC[marker == input$ctrlMarker] != "Ok") "Bad Control Marker"
          else "Ok"
        ) %>%
        group_by(kit) %>% 
        mutate(
          ntc_QC = if (any(cq_QC[sample.type == "ntc"] == "Ok"))
            "Bad NTC" else "Ok",
          kit_QC = ifelse(ntc_QC != "Ok", "Kit Error", "Ok")
        ) %>% 
        group_by(kit, allele) %>% 
        mutate(
          pos_QC = if (all(ctrlMarker_QC[sample.type == "pos"] != "Ok") &
                       xor(all(cq_QC[sample.type == "pos" & 
                                     grepl("ref", sample) & 
                                     marker != input$ctrlMarker] == "Ok") & 
                           all(meanCq_QC[sample.type == "pos" & 
                                         grepl("ref", sample) & 
                                         marker != input$ctrlMarker] == "Ok") &
                           all(deltacq_QC[sample.type == "pos" & 
                                          grepl("ref", sample) & 
                                          marker != input$ctrlMarker] == "Ok"),
                           all(cq_QC[sample.type == "pos" & 
                                     grepl("alt", sample) & 
                                     marker != input$ctrlMarker] == "Ok") & 
                           all(meanCq_QC[sample.type == "pos" & 
                                         grepl("alt", sample) & 
                                         marker != input$ctrlMarker] == "Ok") &
                           all(deltacq_QC[sample.type == "pos" & 
                                          grepl("alt", sample) & 
                                          marker != input$ctrlMarker] == "Ok")))
          "Ok" else "Bad Pos"
          # TODO: check alt = Ok & alt = Ok   or  ref = Ok & ref = Ok
          
        )
      fRes <- res %>% 
        # filter(ctrlMarker_QC == "Ok", kit_QC == "Ok") %>% 
        group_by(kit, marker, sample) %>% 
        mutate(result = paste(allele[cq_QC == "Ok"] %>% unique(), collapse = ""),
               resultZygosity = 
                 sapply(result, 
                        function(res) switch(as.character(str_length(res)),
                                             "0" = "",
                                             "1" = "Homo", "2" = "Hetero", "Error"))
        )
      res <- left_join(res, fRes)
      # testTbl <<- res
      # combine all errors in one column
      res$total_QC <- apply(res %>% 
                              ungroup() %>% 
                              dplyr::select(ends_with("_QC")), 1,
                            function(x) paste(x[x != "Ok"], collapse = ", "))
      res
    })
  })
  
  output$pcrPlateUI <- renderUI({
    req(resultsTbl())
    toLog("Creating pcrPlate")
    pcrPlateInput("pcrPlate", 
                  plateDescription = resultsTbl() %>% 
                    filter(kit %in% input$showKits & marker %in% input$showMarkers) %>% 
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
                  ")
  })
  
  output$ampCurvesUI <- renderUI({
    req(calcResults())
    toLog("Drawing curves")
    renderAmpCurves("ampCurves", 
                    ampCurves = calcResults()$fData,
                    colorBy = "marker")
  })
  
  observeEvent(
    c(input$pcrPlate, input$showMarkers),
    {
      toLog("Updating curves")
      toHideCurves <-
        which(!(calcResults()$dTbl$position %in% input$pcrPlate) |
                !(calcResults()$dTbl$kit %in% input$showKits) |
                !(calcResults()$dTbl$marker %in% input$showMarkers)
              # |
              #   !(calcResults()$dTbl$target.dyeId %in% input$showDyes)
        )
      updateCurves(session,
                   "ampCurves",
                   hideCurves = toHideCurves)
    })
  
  output$globalResultsTbl <- renderDataTable({
    req(resultsTbl())
    toLog("Creating ResultsTbl")
    resultsTbl() %>%
      filter(position %in% input$pcrPlate &
               kit %in% input$showKits &
               marker %in% input$showMarkers) %>% 
      ungroup() %>% 
      dplyr::select(position, marker, kit, sample, sample.type,
                    Cq = cq_f, Mean_Cq = meanCq_f, Delta_Cq = deltaCq_f,
                    result, Zygosity = resultZygosity,
                    kit_QC, total_QC)
  })
})