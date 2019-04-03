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

dataType$public_methods
dataType$set("public", "Preprocess", 
             function(bgRange) {
               # library(chipPCR)
               library(qpcR)
               fpoints <- self$adp$fpoints
               fpoints$fluor <- 
                 chipPCR::CPP(fpoints$cyc, fpoints$fluor,
                     bg.range = bgRange)$y.norm
               self$adp$fpoints <- fpoints
               
               self$endPt <- tail(fpoints$fluor, 1)
               result <- "Ok"
               self$cq <- 
                 tryCatch(
                 efficiency(
                   pcrfit(fpoints[, c("cyc", "fluor")], 
                                cyc = 1, fluo = 2,
                                model = l5),
                   plot = FALSE,
                   type = "cpD2"
                 )$cpD2,
                 error = function(e) {
                   result <<- list(date(), "Calc Cq error")
                   max(fpoints$cyc)
                 })
               result
             },
             overwrite = TRUE)