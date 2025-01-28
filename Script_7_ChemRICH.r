##############
# Script 7: ChemRICH Analysis for Metabolite Set Enrichment
##############
# Clean global environment
rm(list = ls())


# ##############
# # Install Required R Packages (Run once)
# ##############
# if (!require("devtools"))
#     install.packages('devtools', repos="http://cran.rstudio.com/")
# 
#   if (!require("RCurl"))
# 
#     install.packages('RCurl', repos="http://cran.rstudio.com/")
# 
#   if (!require("pacman"))
# 
#     install.packages('pacman', repos="http://cran.rstudio.com/")
# 
#   library(devtools)
# 
#   library(RCurl)
# 
#   library(pacman)
# 
#   if (!requireNamespace("BiocManager", quietly = TRUE))
# 
#     install.packages("BiocManager")
# 
#   pacman::p_load(GGally)
# 
#   pacman::p_load(DT)
# 
#   pacman::p_load(RCurl)
# 
#   pacman::p_load(RJSONIO)
# 
#   pacman::p_load(ape)
# 
#   pacman::p_load(devEMF)
# 
#   pacman::p_load(dynamicTreeCut)
# 
#   pacman::p_load(extrafont)
# 
#   pacman::p_load(ggplot2)
# 
#   pacman::p_load(ggpubr)
# 
#   pacman::p_load(ggrepel)
# 
#   pacman::p_load(grid)
# 
#   pacman::p_load(htmlwidgets)
# 
#   pacman::p_load(igraph)
# 
#   pacman::p_load(magrittr)
# 
#   pacman::p_load(network)
# 
#   pacman::p_load(officer)
# 
#   pacman::p_load(openxlsx)
# 
#   pacman::p_load(phytools)
# 
#   pacman::p_load(plotly)
# 
#   pacman::p_load(plotrix)
# 
#   pacman::p_load(rcdk)
# 
#   pacman::p_load(readxl)
# 
#   pacman::p_load(rvg)
# 
#   pacman::p_load(sna)
# 
#   pacman::p_load(visNetwork)
#   
# # Install Java and rJava manually

 

##############
# Functions
##############
checkSmiles <- function(input = "inputfile.xlsx") {

  ndf <- readxl::read_xlsx(input)

  fps <- lapply(1:nrow(ndf), function(x) {

    rcdk::parse.smiles(ndf$smiles[x])

  })

  charvec <- sapply(fps, nchar)

  paste0("Lines with an incorrect SMILES codes are : ", paste(as.integer(which(charvec==4)), collapse = ","))

}


##############
# Values to Change
##############
wd <- "C:\\Users\\lazab\\Documents\\github\\GF_GCMS_data_analysis"
mesh_prediction_input_filename <- "GF_GCMS_MESH_prediction_input.xlsx"
chemrich_input_filename <- "GF_GCMS_my_batch_ChemRICH_R_input.xlsx"


##############
# Set working directory
##############
output_dir <- paste(wd, "output", sep = "\\")
setwd(output_dir)
wd_output <- getwd()


# ##############
# # Partway through Script 6: Run MESH Predictions (comment out when not in use)
# ##############
# source("https://raw.githubusercontent.com/barupal/ChemRICH/master/predict_mesh_chemical_class.R")
# load.ChemRICH.Packages()
# 
# # Check if there are incorrect SMILES in the input file
# checkSmiles(mesh_prediction_input_filename)
# 
# predict_mesh_classes(inputfile = mesh_prediction_input_filename)
# 
# # Then, manually edit/curate results as needed before returning to Script 6


##############
# Script 7
##############

##############
# ChemRICH: cluster by Tanimoto Score
##############


##############
# ChemRICH: cluster by lipophilicity
##############
source("https://raw.githubusercontent.com/barupal/ChemRICH/master/chemrich_chemical_classes.R")
load.ChemRICH.Packages()

run_chemrich_chemical_classes(inputfile = chemrich_input_filename)


setwd(wd)
# If .Rhistory was created in the working directory, delete it
if (file.exists(".Rhistory")) {
  file.remove(".Rhistory")
}