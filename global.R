# Define globally available objects

# imports -----------------------------------------------------------------

library(shiny)
library(shinycssloaders)
library(shinyBS)
library(shinydashboard)
library(shinythemes)
library(plotly)

library(network)
library(sna)
library(igraph)
library(intergraph)
library(visNetwork)

library(mixOmics)
library(ggmixOmics)

library(tidyverse)
library(cowplot)
library(GGally)
library(ggnetwork)

source('helpers.R')

# data ----
M <- readRDS('data/TCGA model.rds')
corMat <- abs(getCorMat(M))

model1 <- M
model2 <- M

# Get component names
dataNames <- names(M$X)
nEntries <- length(dataNames)
nComp <- unique(M$ncomp)

# Params ----
geneEnrichment <- TRUE
PPIIntegration <- FALSE

quarterWidth <- 3
halfWidth <- 6
tQuarterWidth <- 9
fullWidth <- 12
