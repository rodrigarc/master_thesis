library(Biostrings)
library(stringr)
library(DECIPHER)
library(seqinr)

hc <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/AA sequences/LOR-HC.fasta")
hc_uca <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/AA sequences/LOR-HC-UCA copy.fasta")
lc <- readAAStringSet(
"~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/AA sequences/LOR-LC.fasta")
lc_uca <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/AA sequences/LOR-LC-UCA copy.fasta")

chain_align <- function(HC, HC_uca, LC, LC_uca, x){
heavy_chain <- pairwiseAlignment(HC_uca[x], HC[x], type = "global",
                    gapOpening = 10,
                    gapExtension = 0.5)
light_chain <- pairwiseAlignment(LC_uca[x], LC[x], type = "global",
                    gapOpening = 10,
                    gapExtension = 0.5)
mab <- list(c(heavy_chain, light_chain))
return(mab)
}

x <- chain_align(hc, hc_uca, lc, lc_uca, 1)
x[2]
y <- writePairwiseAlignments(x[1])

seq_hc <- c(alignedPattern(x), alignedSubject(x))
seq_lc <- c(alignedPattern(x[2]), alignedSubject(x[2]))
alignedSubject(x)

seq_hc
seq_lc
x[2]
y <- BrowseSeqs(seq_lc)
pattern(x[1])
pattern(x[2])

global_Align <- pairwiseAlignment(hc_uca[1], hc[1], type = "global",
                                  gapOpening = 10,
                                  gapExtension = 0.5)
pattern(global_Align)
global_Align
global_Align <- list()
for (i in 1:68){
  global_Align[i] <- pairwiseAlignment(hc_uca[i], hc[i],
                                       type = "global",
                                       gapOpening = 10,
                                       gapExtension = 0.5)
}

aligned <- BStringSet()
for (i in 1:68){
test <- BStringSet(c(toString(subject(global_Align[i])), 
                        toString(pattern(global_Align[i]))))
aligned <- append(aligned,test)
}
name <- c(rbind(names(hc), names(hc_uca)))
names(aligned) <- name
score(global_Align)

global_Align[2]

getwd()
setwd("Desktop/")
writeXStringSet(aligned,"out1.txt")

colors <- rainbow(20, v=0.8, start=0.9, end=0.7)
amino <- unique(GENETIC_CODE)
amino <- amino[amino != '*']
amino
m <- match(GENETIC_CODE, unique(GENETIC_CODE)) 
BrowseSeqs(x[1], highlight = 1,
           patterns = c('-',amino),
           colors=c("black", colors))
aligned
class(x[1])
