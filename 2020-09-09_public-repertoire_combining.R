install.packages("ggnewscale")

library(data.table)
library(kableExtra)
library(data.table)
library(Biostrings)
library(dplyr)
library(seqRFLP)
library(openxlsx)

setwd("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire")
#selecting files of clonoquery and reading it
ls <- list.files("public_clonoquery/B2/", full.names = T)
ls <- grep("full", ls, value = T)
tab <- lapply(ls, fread, fill = T, blank.lines.skip = F, 
               stringsAsFactors = F,
               header = T)
#remove empty rows with 0 clones 
#tab_n <- lapply(tab, na.omit)

#combine list into df
tab_c <- rbindlist(tab, fill = T, idcol = T)
ls
write.csv(tab_c,
          "public_clonoquery/B2/B2_combined-clonoquery-full.csv", row.names = F)

#######read csvs and convert to fasta
B1_P0913 <- read.xlsx("../../../Monday Seminar - Sebastian/B1_clonoquery_P0913.xlsx")
B1_P2029 <- read.xlsx("../../../Monday Seminar - Sebastian/B1_clonoquery_P2029.xlsx")
setwd("~/Box Sync/Monday Seminar - Sebastian")
B2_P0254 <- read.xlsx("B2_clonoquery_P0254.xlsx")
IgM_P0254 <- read.xlsx("IgM_clonoquery_P0254.xlsx")

#P0913
B1_P0913 <- na.omit(B1_P0913)
B1_P0913_nt <- B1_P0913 %>% select(name, VDJ_nt)
B1_P0913_aa <- B1_P0913 %>% select(name, VDJ_aa)
dataframe2fas(B1_P0913_aa, "B1_P0913_aa.fasta")
dataframe2fas(B1_P0913_nt, "B1_P0913_nt.fasta")

#P2029
B1_P2029 <- na.omit(B1_P2029)
B1_P2029_nt <- B1_P2029 %>% select(name, VDJ_nt)
B1_P2029_aa <- B1_P2029 %>% select(name, VDJ_aa)

dataframe2fas(B1_P2029_aa, "B1_P2029_aa.fasta")
dataframe2fas(B1_P2029_nt, "B1_P2029_nt.fasta")

#P0254
B2_P0254 <- na.omit(B2_P0254)
B2_P0254_nt <- B2_P0254 %>% select(name, VDJ_nt)
B2_P0254_aa <- B2_P0254 %>% select(name, VDJ_aa)

dataframe2fas(B2_P0254_aa, "B2_P0254_aa.fasta")
dataframe2fas(B2_P0254_nt, "B2_P0254_nt.fasta")


#P0254
IgM_P0254 <- na.omit(IgM_P0254)
IgM_P0254_nt <- IgM_P0254 %>% select(name, VDJ_nt)
IgM_P0254_aa <- IgM_P0254 %>% select(name, VDJ_aa)

dataframe2fas(IgM_P0254_aa, "IgM_P0254_aa.fasta")
dataframe2fas(IgM_P0254_nt, "IgM_P0254_nt.fasta")
