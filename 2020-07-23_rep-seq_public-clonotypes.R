install.packages("R.utils")
library(R.utils)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)
getwd()
#read tab files using fread since the tab files were big
#do it one by one since the computer cannot handle many at the same time 

ls_rs <- list.files("~/Box Sync/public_repertoire/data/2020-09-12/B2_filtered tabs", 
           full.names = T)
ls_rs
#ls_rs <- grep("filtered.tab.gz", ls_rs, value = T)

#tab <- fread(ls_rs[1], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[2], header = T, stringsAsFactors = F)
#tab_3 <- fread(ls_rs[3], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[4], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[5], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[6], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[7], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[8], header = T, stringsAsFactors = F)
#tab <- fread(ls_rs[9], header = T, stringsAsFactors = F)

#create columns based on the ID
df_split <- tab %>% separate(name, sep = "_", c("name"))

#break down alleles to genes in a list - loop
df_rename <- df_split %>% separate(V_gene, into = "V_gene",sep = "\\*|_")
df_rename <- df_rename %>% separate(J_gene, into = "J_gene",sep = "\\*|_")
fwrite(df_rename, 
       file = "~/Box Sync/public_repertoire/data/2020-09-12/B2_processed/E18_B2_edited.tab.gz",
       quote = F, sep = "\t", compress = "gzip")

