library(tidyverse)
library(readxl)    
library(tidyr)
library(reshape)
library(openxlsx)
getwd()
setwd("~/Desktop/rodrigo/V-J Pairing /")

#function to read different worksheets in a excel file and save as a list
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # pass tibble = TRUE if you think tibble is better
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#reading data > although I used the data from the code below for the VJ pairing
ls <- read.xlsx("~/Desktop/rodrigo/V-J Pairing /List of V and J genes.xlsx")
ls <- ls %>% dplyr::rename(V_gene = Unique.Genes, J_gene = J.Genes)
vj <- read.xlsx("~/Desktop/rodrigo/V-J Pairing /VJ frequencies_TR1.xlsx")
vj <- vj %>% dplyr::rename(J_gene = J.gene)

#using the function created above and setting proper column names
data_list <- read_excel_allsheets("~/Desktop/rodrigo/V-J Pairing /VJ frequencies_TR1_individual.xlsx")
new_col_name <- colnames(vj)
data_list <- lapply(data_list, setNames, nm = new_col_name)

#creating every possible unique combination of V, J genes, and specificity
ls$Specificity <- c("DP","PreF")
ls_comb <- expand.grid(ls$V_gene, na.omit(ls$J_gene), na.omit(ls$Specificity))
ls_comb <- ls_comb %>% distinct(.keep_all = TRUE) 
colnames(ls_comb) <- colnames(ls)
names(data_list)


#finally, this one worked to combine data, ignore PostF, 
#calculate sum, and add in a list
x <- data_list[[1]]
merged_x <- merge(x, ls_comb, all = T) %>%
  filter(Specificity != "PostF") %>%
  mutate(ID = unique(na.omit(x$ID)), Size = ifelse(is.na(Size),0,Size)) %>%
  group_by(Specificity,V_gene,J_gene) %>%
  mutate(sum = sum(Size)) %>%
  select(ID, Specificity, V_gene, J_gene, sum) %>%
  ungroup() %>%
  distinct(.keep_all = TRUE) %>%
  group_by(ID, Specificity) %>% 
  group_split() 

#making a loop with the code above
df_list <- list()
for (i in names(data_list)) {
  df_list[[i]] <- merge(data_list[[i]], ls_comb, all = T) %>%
    filter(Specificity != "PostF") %>%
    mutate(ID = unique(na.omit(data_list[[i]][["ID"]])), Size = ifelse(is.na(Size),0,Size)) %>%
    group_by(Specificity,V_gene,J_gene) %>%
    mutate(sum = sum(Size)) %>%
    select(ID, Specificity, V_gene, J_gene, sum) %>%
    ungroup() %>%
    distinct(.keep_all = TRUE) %>%
    group_by(Specificity) %>% 
    group_split() 
}
#unlist 
df_unlist <- unlist(df_list, recursive = FALSE)


#creating pivot tables
df_pivot <- list ()

for (i in names(df_unlist)){
  df_pivot[[i]] <- df_unlist[[i]] %>% 
    cast(V_gene ~ J_gene)
}
names(df_pivot) <- c("E24_DP","E24_PreF",
                     "E23_DP","E23_PreF",
                     "E21_DP","E21_PreF",
                     "E18_DP","E18_PreF",
                     "E17_DP","E17_PreF",
                     "E16_DP","E16_PreF",
                     "E14_DP","E14_PreF",
                     "E12_DP","E12_PreF",
                     "E11_DP","E11_PreF")

#saving as a worksheet
write.xlsx(df_pivot,"V_J.xlsx")


#do the same, but for total
df_list_total <- list()
for (i in names(data_list)) {
  df_list_total[[i]] <- merge(data_list[[i]], ls_comb, all = T) %>%
    mutate(ID = unique(na.omit(data_list[[i]][["ID"]])), Size = ifelse(is.na(Size),0,Size)) %>%
    group_by(V_gene,J_gene) %>%
    mutate(sum = sum(Size)) %>%
    select(ID, V_gene, J_gene, sum) %>%
    ungroup() %>%
    distinct(.keep_all = TRUE) %>%
    group_split() 
}
#unlist 
df_unlist_total <- unlist(df_list_total, recursive = FALSE)

#creating pivot tables
df_pivot_total <- list ()
for (i in names(df_unlist_total)){
  df_pivot_total[[i]] <- df_unlist_total[[i]] %>% 
    cast(V_gene ~ J_gene)
}


#saving all (total, pref, DP) in a worksheet
final_list <- c(df_pivot_total, df_pivot)
#saving as a worksheet
write.xlsx(final_list,"V_J_pairing.xlsx")

#calculate the percentages for each pivot table
#bind them and create a new xlsx file 
df_pivot_percent <- lapply(df_pivot, function(x) prop.table((as.matrix(x)))*100)
df_pivot_total_percent <-  lapply(df_pivot_total, function(x) prop.table((as.matrix(x)))*100) 
final_list_percent <- c(df_pivot_total_percent, df_pivot_percent)
write.xlsx(final_list_percent,"V_J_percentage.xlsx")

