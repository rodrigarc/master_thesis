library(Biostrings)


ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined", 
                 full.names = T)
s1 <- readDNAStringSet(ls[1])
s2 <- readDNAStringSet(ls[2])
s3 <- readDNAStringSet(ls[3])
s4 <- readDNAStringSet(ls[4])
s5 <- readDNAStringSet(ls[5])
s6 <- readDNAStringSet(ls[6])
s7 <- readDNAStringSet(ls[7])
s8 <- readDNAStringSet(ls[8])
s9 <- readDNAStringSet(ls[9])

fs <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9)
fs <- translate(fs)
nm <- names(fs)
nm_u <- unique(nm)
fs_unique <- fs[nm_u]
fs_unique

writeXStringSet(fs_unique, "~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined/v_aa.fasta")



#add "-" as many times as it appears
library(xlsx)
library(stringr)
library(dplyr)

dsh <- read.xlsx("~/Downloads/LOR HCDR3 dashes.xlsx", sheetIndex = 1)
setwd("~/Desktop/")

dsh <- dsh %>% mutate("-" = strrep("-",dsh$X12))


write.xlsx(dsh, "LOR HCDR3 dashes.xlsx")








