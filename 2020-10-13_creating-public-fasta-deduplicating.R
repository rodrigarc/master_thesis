library(data.table)
library(seqRFLP)
library(Biostrings)




######selecting public clones #######
df_comb_selection <- df_comb_jit_with_LC %>% filter(number_animals_combined >= 5)
p_clone_selection <- unique(df_comb_selection$clone)
p_clone_sequences <- combined_full[combined_full$clone %in% p_clone_selection]

p_nt <- p_clone_sequences %>% select(name, VDJ_nt, clone) 
p_nt$VDJ_nt <- gsub("\\*","", p_nt$VDJ_nt)
w <- "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/not_dedup/"

for (i in unique(p_nt$clone)){
dataframe2fas(p_nt[p_nt$clone == i, c("name","VDJ_nt")], paste0( w, i,"_public_clones_5animals_nt.fasta"))
}

p_aa <- p_clone_sequences %>% select(name, VDJ_aa, clone)
p_aa$VDJ_aa <- gsub("\\*","", p_aa$VDJ_aa)

for (i in unique(p_aa$clone)){
  p <- p_aa[p_aa$clone == i, c("name","VDJ_aa")]
  dataframe2fas(p, paste0( w, i,"_public_clones_5animals_aa.fasta"))
}
View(table(p_aa$VDJ_aa[p_aa$clone == "P2029"]))

ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/not_dedup/",full.names = T)
ls_aa <- grep("aa", ls, value = T)
ls_nt <- grep("nt", ls, value = T)


fasta_aa <- sapply(ls_aa,readBStringSet,simplify = FALSE,USE.NAMES = TRUE)
fasta_nt <- sapply(ls_nt,readBStringSet,simplify = FALSE,USE.NAMES = TRUE)

fasta_aa_u <- lapply(fasta_aa, unique)
names(fasta_aa_u) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")
fasta_nt_u <- lapply(fasta_nt, unique)
names(fasta_nt_u) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")

z <- "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/dedup/"

for (i in names(fasta_aa_u)){
  writeXStringSet(fasta_aa_u[[i]], paste0(z, i, "_dedup_public_clones_5animals_aa.fasta"))
}

for (i in names(fasta_nt_u)){
  writeXStringSet(fasta_nt_u[[i]], paste0(z, i, "_dedup_public_clones_5animals_nt.fasta"))
}




