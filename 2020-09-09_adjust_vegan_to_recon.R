#code to adjust vegan file to recon.py

library(dplyr)
getwd()
setwd("~/Desktop/Recon/2020-09-09")

#delete NA
dp[is.na(dp)] <- 0
pref[is.na(pref)] <- 0
total[is.na(total)] <- 0

#transpose
dp_t <- t(dp)
pref_t <- t(pref)
total_t <- t(total)

#split matrix columns into new lists
dp_t_split <- asplit(dp_t, 2)
pref_t_split <- asplit(pref_t, 2)
total_t_split <- asplit(total_t, 2)


#remove 0s from tables to recalculate (may skip this step if do not want to delete 0)
dp_t_drop <- lapply(dp_t_split, function(x) {x[x!=0]})
pref_t_drop <- lapply(pref_t_split, function(x) {x[x!=0]})
total_t_drop <- lapply(total_t_split, function(x) {x[x!=0]})

##count how many times each clone without 0s
dp_t_table_drop <- lapply(dp_t_drop, table)
pref_t_table_drop <- lapply(pref_t_drop, table)
total_t_table_drop <- lapply(total_t_drop, table)

#count how many times each clone with 0s
dp_t_table <- lapply(dp_t_split, table)
pref_t_table <- lapply(pref_t_split, table)
total_t_table <- lapply(total_t_split, table)


#export as tab files for dp
for(i in 1:length(dp_t_table_drop)) {
  write.table(dp_t_table_drop[[i]], 
              file = paste0(i,"_dp", ".txt"), 
              row.names = F,
              col.names = F, 
              quote = F,
              sep = "\t")
}

#export as tab files for pref
for(i in 1:length(pref_t_table_drop)) {
  write.table(pref_t_table_drop[[i]], 
              file = paste0(i,"_pref", ".txt"), 
              row.names = F,
              col.names = F, 
              quote = F,
              sep = "\t")
}

#export as tab files for total
for(i in 1:length(total_t_table_drop)) {
  write.table(total_t_table_drop[[i]], 
              file = paste0(i,"_total", ".txt"), 
              row.names = F,
              col.names = F, 
              quote = F,
              sep = "\t")
}

#BASH CODE FOR RUNNING IN TERMINAL ALL OF THE PYTHON SCRIPT IN ALL FILES IN A FOLDER
for f in *.txt; do
 python2 recon_v2.5.py -R -t 30 -c -o "${f%.txt}_fitfile.txt" "$f" 
 done
 
#BASH CODE FOR r SUBSAMPLE
 for f in *.txt; 
 do python2 ../../recon_v2.5.py -r -o "${f%.txt}_resample.txt" "$f" ; 
 done

#bash code for plotting
 for f in *fitfile.txt;
 do python2 ../../recon_v2.5.py -x -o "${f%.txt}_resample.txt" -b error_bar_parameters.txt "$f" ;
 done
