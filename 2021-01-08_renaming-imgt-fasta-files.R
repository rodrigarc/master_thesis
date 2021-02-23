library(DECIPHER)
library(stringr)
library(Biostrings)

##### editing IMGT file names
lc_k <- readDNAStringSet("~/Box Sync/rsv-sc/results/2020-12-17/IMGThs_LC/IGKV.fasta") #reading fasta
lc_k <- RemoveGaps(lc_k) ## removing gaps, "." and "-"
lc_k_name <- names(lc_k)
lc_k_name
rename_k <- str_split(lc_k_name,"\\|", n = 3) ## splitting names based on "|"
rename_k <- sapply(rename_k,FUN = function(x) x[[2]]) ## take only the allele names, the second match
rename_k
names(lc_k) <- rename_k ## renaming
lc_k

##### editing IMGT file names
lc_l <- readDNAStringSet("~/Box Sync/rsv-sc/results/2020-12-17/IMGThs_LC/IGLV.fasta") #reading fasta
lc_l <- RemoveGaps(lc_l) ## removing gaps, "." and "-"
lc_l_name <- names(lc_l)
lc_l_name
rename_l <- str_split(lc_l_name,"\\|", n = 3) ## splitting names based on "|"
rename_l <- sapply(rename_l, FUN = function(x) x[[2]]) ## take only the allele names, the second match
rename_l
names(lc_l) <- rename_l ## renaming
lc_l


#### combining kappa and lambda
light_chain <- c(lc_k, lc_l)
light_chain

writeXStringSet(light_chain, "~/Box Sync/rsv-sc/results/2020-12-17/IMGThs_LC/V.fasta")
