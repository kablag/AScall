toLog <- function(logMessage, timeStamp =  date()) {
  cat(sprintf("%s: %s \n", timeStamp, logMessage))
}

addMarkersKits <- function(descrTbl) {
  descrTbl %>% 
    mutate(marker = str_split(target, "_", simplify = TRUE)[, 1],
           allele = str_split(target, "_", simplify = TRUE)[, 2]) %>% 
    group_by(position) %>% 
    mutate(kit = paste(marker, collapse = ", ")) %>% 
    ungroup() %>% 
    mutate(kit = kit %>% 
             as.factor() %>% 
             as.numeric() %>% 
             sprintf("kit_%02i", .))
}

dataType$set("public", "rawAdp",
             NA,
             overwrite = TRUE)
dataType$set("public", "rawCq",
             NA,
             overwrite = TRUE)

dataType$set("public", "Init",
             function() {
               if (!is.environment(self$rawAdp)) {
                 # first preprocess init
                 self$rawAdp <- self$adp$copy()
                 self$rawCq <- self$cq
               } else {
                 self$adp <- self$rawAdp$copy()
                 self$cq <- self$rawCq
               }
             },
             overwrite = TRUE)

dataType$set("public", "Preprocess", 
             function(bgRange, modelType, cqMethod, thrCq) {
               # library(chipPCR)
               library(qpcR)
               fpoints <- self$adp$fpoints
               fpoints$fluor <- 
                 chipPCR::CPP(fpoints$cyc, fpoints$fluor,
                              trans = TRUE,
                              bg.range = bgRange)$y.norm
               self$adp$fpoints <- fpoints
               result <- "Ok"
               if (thrCq > max(fpoints$fluor)) 
                 thrCq <- max(fpoints$fluor)
               tryCatch({
                 ampModel <- pcrfit(fpoints[, c("cyc", "fluor")], 
                                    cyc = 1, fluo = 2,
                                    model = get(modelType))
                 fpoints$fluor <- predict(ampModel)
                 self$adp$fpoints <- fpoints
                 self$cq <-
                   efficiency(
                     ampModel,
                     plot = FALSE,
                     type = if (cqMethod == "Threshold") "cpD2" else cqMethod,
                     thresh = if (cqMethod == "Threshold") thrCq else NULL
                   )[[
                     {
                       switch(cqMethod,
                              "cpD2" = "cpD2",
                              "cpD1" = "cpD1",
                              "maxE" = "cpE",
                              "expR" = "cpR",
                              "CQ" = "cpCQ",
                              "Cy0" = "Cy0",
                              "Threshold" = "cpT"
                       )
                     }
                     ]]
               },
               error = function(e) {
                 result <<- list(date(), "Calc Cq error")
                 max(fpoints$cyc)
               })
               result
             },
             overwrite = TRUE)


dataType$set("public", "InitEndPt", 
             function() {
               self$endPt <- tail(self$adp$fpoints$fluor, 1)
             },
             overwrite = TRUE)

runType$set("public", "Preprocess",
            function(bgRange, modelType, cqMethod, thrCq) {
              self$react <- foreach(reactN = self$react) %dopar% {
                for (dat in reactN$data)
                  dat$Preprocess(bgRange, modelType, cqMethod, thrCq)
                reactN
              } 
            },
            overwrite = TRUE)

runType$set("public", "Init",
            function(bgRange, modelType, cqMethod, thrCq) {
              self$react <- foreach(reactN = self$react) %dopar% {
                for (dat in reactN$data)
                  dat$Init()
                reactN
              } 
            },
            overwrite = TRUE)

