
library(Biostrings)
library(msa)
library(ggtree)
library(ape)
library(dplyr)


fs <- readDNAStringSet("~/Box Sync/RSV NGS/v4_Analysis/LOR21_lineage/E16_2_allseqs_nt.fasta")
fs_duplicated <- fs[duplicated(fs)]
fs_unique <- unique(fs)
fs_duplicated
fs_unique

#writeXStringSet(fs_unique, "~/Box Sync/RSV NGS/v4_Analysis/LOR21_lineage/E16_2_allseqs_nt_unique.fasta")

duplicates <-  match(fs_duplicated, fs_unique)
counts <- data.frame(table(duplicates))
counts
n <- as.numeric(counts$duplicates)
counts$names <- names(fs_unique[n])


counts_s <- counts %>% select(names, Freq) 
#write.csv(counts_s,"~/Box Sync/RSV NGS/v4_Analysis/LOR21_lineage/duplicates.csv",row.names = F)

x <- msaMuscle(fs, verbose = T)
aligned_bin <- as.DNAbin(x)
d <- dist.dna(aligned_bin)
tree <- nj(d)
g <- ggtree(tree, layout = "rectangular", MAX_COUNT = 1)
g + geom_tiplab(size = 0.8) 
g1 <- g %<+% counts_s
g1 + geom_tiplab(size = 0.8) +  geom_tippoint(aes(size=Freq))
