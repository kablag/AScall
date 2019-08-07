#list of packages required
list.of.packages <- c("chipPCR", "qpcR", "DT", "doParallel",
                      "tidyverse", "plotly", "RColorBrewer",
                      "shiny", "RDML", "shinyMolBio",
                      "shinyWidgets", "openxlsx")

#checking missing packages from list
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

#install missing ones
if (length(new.packages)) {
  cat(sprintf("Installing missing packages: %s\n", paste(new.packages, collapse = ", ")))
  install.packages(new.packages, dependencies = TRUE)
}

library(RDML)
library(tidyverse)
library(plotly)
library(doParallel)
library(RColorBrewer)
library(openxlsx)

library(shiny)
library(shinyMolBio)
library(shinyWidgets)

nCores <- detectCores()
cl <- makeCluster(nCores)
registerDoParallel(cl)