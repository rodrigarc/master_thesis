#Cross-reference light chains from public and all
install.packages("sjmisc")
library(dplyr)
library(openxlsx)
library(stringr)
library(sjmisc)

#reading files 
LC_all <- read.xlsx("~/Box Sync/RSV NGS/v4_Analysis/LCs/200806_LC_all.xlsx")
LC_public <- read.xlsx("~/Box Sync/RSV NGS/v4_Analysis/LCs/200828_LC_public.xlsx") 
LC_B2 <- read.csv(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/public-tables/2020-08-16_public-table-combined.csv")
LC_B2_only <-  read.csv("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/public-tables/2020-09-13_public_table_B2_summary.csv")


#selecting just clones in B2
LC_B2 <- LC_B2[LC_B2$B2 == 1,]

#extracting just the well names and ignore chain and specificity
LC_B2$well <- sapply(strsplit(LC_B2$well, "_"), function(x) paste(x[c(-5,-4)], collapse = "_"))
LC_public <- LC_public %>% rename(well = name)

#merging dfs
LC_B2_public <- merge(LC_B2, LC_public, by = "well")

###################################################################################
#greping the sequences names from LC_public within B2
#tried package sjmisc, the best vectorized pattern matching I found
#but at end I just used loops and grep...

toyData <- LC_B2_only  
vars <- LC_all$well     
# create a vector to collect the output of each call    

toyColIndexList <- vector(length = length(vars), mode = "list")    

# grep each element in turn     

for (i in seq_along(vars)) {      
  toyColIndexList[[i]] <- grep(pattern = vars[i], x = toyData$sequence)     
}      

# combine all of the answers     
toyColIndex <- unlist(toyColIndexList)   

# remove duplicated columns if present    
toyColIndex <- toyColIndex[!duplicated(toyColIndex)]     

# select the elements we want    

B2_public <- toyData[toyColIndex,]  
B2_all <- toyData[toyColIndex,] 

write.xlsx(B2_public,
           "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/public-tables/2020-09-15_LC_B2-public.xlsx")
write.xlsx(B2_all,
           "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/public-tables/2020-09-15_LC_B2-all.xlsx")








