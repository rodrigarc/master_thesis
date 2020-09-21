library(tidyverse)
library(readxl)    
library(tidyr)
library(reshape)
library(openxlsx)
library(stringr)
getwd()
setwd("~/Desktop/rodrigo/V-J Pairing /")

#reading files, putting in a list, deleting extrafiles, setting proper column names
load("~/Desktop/IgM_v4_clean.RData")
data_list <- list(E11.IgM, E12.IgM, E14.IgM, E16.IgM, 
                  E17.IgM, E18.IgM, E21.IgM, E23.IgM, E24.IgM)
names(data_list) <- c("E11.IgM", "E12.IgM", "E14.IgM", "E16.IgM", 
                      "E17.IgM", "E18.IgM", "E21.IgM", "E23.IgM", "E24.IgM")
rm(E11.IgM, E12.IgM, E14.IgM, E16.IgM, 
   E17.IgM, E18.IgM, E21.IgM, E23.IgM, E24.IgM)
data_list_select <- lapply(data_list, select, Animal, V_gene, J_gene, count)
data_list_select <- lapply(data_list_select, dplyr::rename, ID = Animal, Size = count)


#reading data > although I used the data from the code below for the VJ pairing
ls <- read.xlsx("~/Desktop/rodrigo/V-J Pairing /List of V and J genes.xlsx")
ls <- ls %>% dplyr::rename(V_gene = Unique.Genes, J_gene = J.Genes)
ls_comb <- expand.grid(ls$V_gene, na.omit(ls$J_gene))
ls_comb <- ls_comb %>% distinct(.keep_all = TRUE) 
colnames(ls_comb) <- colnames(ls)

#break down alleles to genes in a list - loop
data_list_select <- lapply(data_list_select, 
                           separate, V_gene, into = "V_gene",sep = "\\*|_")
data_list_select <- lapply(data_list_select, 
                           separate, J_gene, into = "J_gene",sep = "\\*|_")
data_list_select <- lapply(data_list_select, function (x) mutate(x, J_gene = 
      ifelse(x$J_gene == "IGHJ2", "IGHJ2-1", x$J_gene)))


#calculate sum, and add in a list
#making a loop with the code for all lists
df_list <- list()
for (i in names(data_list_select)) {
  df_list[[i]] <- merge(data_list_select[[i]], ls_comb, all = T) %>%
    mutate(Size = ifelse(is.na(Size),0,Size)) %>%
    group_by(V_gene,J_gene) %>%
    mutate(sum = sum(Size)) %>%
    select(ID, V_gene, J_gene, sum) %>%
    ungroup() %>%
    distinct(.keep_all = TRUE) 
}

#creating pivot tables
df_pivot <- list ()
for (i in names(df_list)){
  df_pivot[[i]] <- df_list[[i]] %>% 
    cast(V_gene ~ J_gene)
}
names(df_pivot) <- c("E11.IgM", "E12.IgM", "E14.IgM", "E16.IgM", 
                     "E17.IgM", "E18.IgM", "E21.IgM", "E23.IgM", "E24.IgM")

#saving as a worksheet
write.xlsx(df_pivot,"V_J_IgM.xlsx")

#calculate the percentages for each pivot table
#bind them and create a new xlsx file 
df_pivot_percent <- lapply(df_pivot, function(x) prop.table((as.matrix(x)))*100)
write.xlsx(df_pivot_percent,"V_J_IgM_percentage.xlsx")
