###### code to count number of alleles discovered by our improved database
library(Biostrings)
comb_v <- readDNAStringSet(filepath = "~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined/combined.fasta")
ws <- readDNAStringSet(filepath = "~/Box Sync/RSV NGS/Database comparison/Cirelli/database/V.fasta")
comb_v
ws

dif <- setdiff(comb_v, ws)
un <- union(comb_v, ws)

length(names(dif))
length(grep("_S",names(dif)))

