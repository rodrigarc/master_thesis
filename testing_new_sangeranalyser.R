install_github("roblanf/sangeranalyseR", ref = "develop")
library(sangeranalyseR)
library(shiny)
library(tidyverse)
getwd()
setwd("~/Desktop/rodrigo/sanger_HC")

sangerReadF <-  SangerAlignment("sanger_test/")

                           
sangerreadF <- SangerRead(inputSource           = "ABIF",
           readFeature           = "Forward Read",
           readFileName          = "sanger_test/A8-3_IgG_Inner_GTTCAGGGAAGTAGTCCTTGAC.ab1",
           geneticCode           = GENETIC_CODE,
           TrimmingMethod        = "M1",
           M1TrimmingCutoff      = 0.0001,
           M2CutoffQualityScore  = NULL,
           M2SlidingWindowSize   = NULL,
           baseNumPerRow         = 100,
           heightPerRow          = 200,
           signalRatioCutoff     = 0.33,
           showTrimmed           = TRUE)

sangerContig <- SangerContig(inputSource = "ABIF",parentDirectory= "sanger_test")
qualityBasePlot(sangerreadF)
sangeranalyseR::generateReport(sangerreadF)
