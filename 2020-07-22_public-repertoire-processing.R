#public clonotypes
library(dplyr)
library(tidyr)
library(openxlsx)


#read all filtered files in a list
ls <- list.files("~/Desktop/rodrigo/public-clonotype/sc_v4 filtered files", 
                 full.names =  TRUE)
ls_tab <- lapply(ls, read.table, header = T, sep = "\t", fill = TRUE)

#combine all of them in one single df
df_all <- do.call(rbind, ls_tab)

#create columns based on the ID
df_split <- df_all %>% separate(name,sep = "_", c("name"))

#reading data > although I used the data from the code below for the VJ pairing
ls_g <- read.xlsx("~/Desktop/rodrigo/V-J Pairing /List of V and J genes.xlsx")
ls_g <- ls_g %>% dplyr::rename(V_gene = Unique.Genes, J_gene = J.Genes)
ls_gcomb <- expand.grid(ls_g$V_gene, na.omit(ls_g$J_gene))
ls_gcomb <- ls_gcomb %>% distinct(.keep_all = TRUE) 
colnames(ls_gcomb) <- colnames(ls_g)

#break down alleles to genes in a list - loop
df_rename <- df_split %>% separate(V_gene, into = "V_gene",sep = "\\*|_")
df_rename <- df_rename %>% separate(J_gene, into = "J_gene",sep = "\\*|_")

#just in case it is needed - but I did not use it > rename J gene 2 
data_list_select <- lapply(data_list_select, function (x) mutate(x, J_gene = 
                                                                   ifelse(x$J_gene == "IGHJ2", "IGHJ2-1", x$J_gene)))

write.table(df_rename, 
            file ="~/Desktop/rodrigo/public-clonotype/filtered.tab", 
            sep = "\t",
            row.names = F,
            col.names = T,
            quote = F)
