library(Biostrings)

ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined", full.names = T)

ls <- ls[-c(1,11)]
dna <- lapply(ls, readDNAStringSet)
dna_c <- lapply(dna, c)
dna_c <- do.call(c,)

id <- c("E11","E12","E14","E16","E17","E18","E21","E23","E24")

ls <- list()
for(i in id) {
 ls[[i]] <- list.files(paste0("~/Box Sync/RSV NGS/Database comparison/individualized/",i,"/final"), 
                       full.names = T)
 ls[[i]] <- grep("filtered", ls[[i]], value = T)
  tabs_ind <- lapply(ls, fread)
}

tabs_comb <- rbindlist(tabs_ind, idcol = T) 
tabs_comb_sel <- tabs_comb %>% select(.id, name, V_SHM)
fwrite(tabs_comb, "~/Box Sync/RSV NGS/Database comparison/filtered_tabs/individualized.tab",
             quote = F)


