library(Biostrings)
library(openxlsx)
library(msa)

combined_full %>% 
  filter(clone %in% c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")) %>% 
  group_by(clone, V_gene, J_gene) %>% summarise() %>% View()

table <- read.xlsx("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-15/public-clones_genes-alleles_query.xlsx", sheet = 1)

ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined", full.names = T)
ls <- grep("aa", ls, value = T, invert = T)
ls_bs <- lapply(ls, readDNAStringSet)

all <- do.call(c,ls_bs)

selection <- all[table$V_alleles[table$clone != "P2029"]] 

##### I WAS GOING TO DO IT IN A LOOP BUT SINCE SEBASTIAN WANTS ITS ULTRAFAST, LET'S GO SIMPLE AND DUMB ####
#selection_ls <- list()
#for(i in unique(table$clone)){
#  selection_ls[[i]] <- selection[grep("IGHV5-15\\*01", names(selection), value = T, fixed = F)]
#}
  
selection_ls <- list()
selection_ls[["P0049"]] <- selection[grep("IGHV5-15\\*01", names(selection), value = T, fixed = F)]
selection_ls[["P0108"]] <- selection[grep("IGHV5-157\\*0", names(selection), value = T, fixed = F)]
selection_ls[["P0254"]] <- selection[grep("IGHV4-NGC30", names(selection), value = T, fixed = F)]
selection_ls[["P0733"]] <- selection[grep("IGHV2-118\\*01", names(selection), value = T, fixed = F)]
selection_ls[["P0741"]] <- selection[grep("IGHV2-173\\*01", names(selection), value = T, fixed = F)]
selection_ls[["P0867"]] <- selection[grep("IGHV3-172\\*01", names(selection), value = T, fixed = F)]
selection_ls[["P0913"]] <- selection[grep("IGHV3-50\\*0", names(selection), value = T, fixed = F)]
selection_ls[["P1039"]] <- selection[grep("IGHV3-94\\*01", names(selection), value = T, fixed = F)]
selection_ls[["P1123"]] <- selection[grep("IGHV4-117\\*01", names(selection), value = T, fixed = F)]
selection_ls

selection_aln <- lapply(selection_ls, msaClustalOmega)

###### this thing downhere worked to sink out the consensus sequence, thanks god #####
printSplitString <- function(x, width=getOption("width") - 1)
{
  starts <- seq(from=1, to=nchar(x), by=width)
  for (i in 1:length(starts))
    cat(substr(x, starts[i], starts[i] + width - 1), "\n")
}
names(selection_aln)
for(i in names(selection_aln)){
  sink(paste0("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-15/UCA/", i,".fasta"))
  cat(paste0('>', i, "_", "consensus"))
  cat('\n')
  the.seq <- toString(printSplitString(msaConsensusSequence(selection_aln[[i]])))
  cat(the.seq)
  cat('\n')
  } 
sink(NULL)
########## FAILLING TO USE OTHER FUNCTIONS
    detail(selection_aln[[1]])
.alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  sink(NULL)
}
for(i in names(selection_aln)){
  writeXStringSet(selection_aln[[i]], 
                  filepath = paste0(
                    "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-15/UCA/", i,".fasta"))
}


>IGHJ5-5*01
ACAACTCATTGGATGTCTGGGGCCGGGGAGTTCTGGTCACCGTCTCCTCAG
>IGHJ5-4*02
ACAACTGGTTCGATGTCTGGGGCCCGGGAGTCCTGGTCACCGTCTCCTCAG

x <- DNAStringSet(c("IGHJ5-5*01" = "ACAACTCATTGGATGTCTGGGGCCGGGGAGTTCTGGTCACCGTCTCCTCAG",
             "IGHJ5-4*02" = "ACAACTGGTTCGATGTCTGGGGCCCGGGAGTCCTGGTCACCGTCTCCTCAG"))
y <- pairwiseAlignment(x[1],x[2])

consensusString(y)




              