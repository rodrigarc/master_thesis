#essential devtools for ubuntu 
#this code below should be run in the terminal
sudo apt install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
######

install.packages("devtools")
devtools::install_github("hadley/devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")
BiocManager::install("sangerseqR")

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("devtools")
library("DECIPHER")
library("Biostrings")
library("sangerseqR")
library("dplyr")
library("ggplot2")
library("tidyverse")
library("gridExtra")

sf <- summarise.abi.folder("sanger_HC", trim.cutoff = 0.01)
