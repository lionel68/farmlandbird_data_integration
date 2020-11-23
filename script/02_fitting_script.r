#################
# Example script to reproduce 
# analysis from:
# Hertzog et al Model-based data integration
# problems can be reported by mailing:
# lionel.hertzog@thuenen.de
#############

# note that due to data sharing policy
# only artificial data are provided 
# in the repo
# please contact the main author
# for request of sharing the actual data

# set the working directory to the 
# local file path of the repo
setwd("")

# load the libraries
library(tidyverse)
library(rtrim)
library(greta)
library(greta.checks)
library(greta.multivariate)

## load functions
source("script/01_function_scripts.r")

## where to save the model outputs
path_out <- "model_output/"

## load the example data
mhb <- read.csv("data/example_mhb.csv", stringsAsFactors = FALSE)

ornitho_list <- list_all <- read.csv("data/example_ornitho_list.csv", stringsAsFactors = FALSE)

ornitho_unstr <- unstr_all <- read.csv("data/example_ornitho_unstr.csv", stringsAsFactors = FALSE)

## loop across the species
for(species in unique(mhb$Artname)){
  # for saving the model objects check that the species has a folder existing
  if(!dir.exists(paste0(path_out, species))){
    dir.create(paste0(path_out, species))
  }
  print(paste0("Starting species: ", species))
  # fit the models
  fit_model1(species, path = path_out)
  fit_model2(species, path = path_out, one_by_one = TRUE)
  fit_model3(species, path = path_out, one_by_one = TRUE)
  fit_model4(species, path = path_out, one_by_one = TRUE)
  fit_model5(species, path = path_out, one_by_one = TRUE)
  fit_model6(species, path = path_out, one_by_one = TRUE)
}

