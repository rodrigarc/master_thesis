install.packages("stringdist")
install.packages("qgraph")
install.packages("pvclust")
install.packages("MASS")
install.packages("colordistance")
library(dplyr)
library(stringdist)
library(Biostrings)
library(qgraph)
library(pvclust)
library(MASS)
library(ggplot2)
library(colordistance)
library(gplots)
library(ggmsa)


v4_V_combined <- readDNAStringSet("/Users/karlor/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined/combined.fasta")
u_v4 <- unique(v4_V_combined)
length(u_v4)

v4_V_combined_sample <- v4_V_combined[sample(1:length(v4_V_combined), 200)]

dat <- data.frame()
for(i in c(10,20,50,100,200)){
  start_time <- Sys.time()
  test <- stringDist(v4_V_combined[sample(1:length(v4_V_combined), i)], method = "levenshtein", ignoreCase = TRUE, 
           diag = FALSE,
           upper = FALSE,
           type = "overlap")
  test_mi <- 1/test
  test_mi_graph <- qgraph(test_mi, layout='spring', vsize=2)
  end_time <- Sys.time()
  df_time <- data.frame(seconds = end_time - start_time, sequences = i)
  dat <- rbind(dat, df_time)
}

#plotting time scalability of sequences
plot(dat)

plot(test_mi_graph)


## testing 200 sequences
test <- stringDist(v4_V_combined[sample(1:length(v4_V_combined), 700)], method = "levenshtein", ignoreCase = TRUE, 
                   diag = FALSE,
                   upper = FALSE,
                   type = "overlap")
test_mi <- 1/test
test_mi_graph <- qgraph(test_mi, layout='spring', vsize=2)

par(cex=1, mar=c(3, 3, 3, 3))
plot(hist(test))



# Ward Hierarchical Clustering
test_200 <- stringDist(v4_V_combined[sample(1:length(v4_V_combined), 200)], method = "levenshtein", ignoreCase = TRUE, 
                   diag = FALSE,
                   upper = FALSE,
                   type = "overlap")
test_mi <- 1/test_200
test_mi_graph <- qgraph(test_mi, layout='spring', vsize=2)

d <- test # distance matrix
fit <- hclust(d, method="average")
plot(fit, cex.main=0.1) # display dendogram
groups <- cutree(fit, k=5) # cut tree into 5 clusters
?cutree
# draw dendogram with red borders around the 5 clusters
rect.hclust(fit, k=5, border="red")

pdf("~/Box Sync/rsv-sc/results/2021-02-12/2021-02-12_hierachical-clustering-3.pdf", 
      width = 12, height =  9)
par(cex=0.1, mar=c(25, 25, 25, 25))
plot(fit, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(xlab="", ylab="", main="Levenstain Distance - Average")
axis(2)
groups <- cutree(fit, h=)
rect.hclust(fit, h=10, border="red")
max(groups)
dev.off()
#create a plot for different number of cuts and count number of clusters
#do a plot for that



#do the same for aminoacid level > fasta file for aminoacids in the mAb=QC combined
#v4_germline > summary germline > fasta files there 

##cluster together aa and nt 

## check the number of sequences we had it - around 80-90 clusters, 60 for aa, redo aa for
#pairwisecomparisons and everything


### MASS

mds <- isoMDS((test_mi))



### running with unique
test <- stringDist(u_v4, method = "levenshtein", ignoreCase = TRUE, 
                   diag = FALSE,
                   upper = FALSE,
                   type = "overlap")
test_mi <- 1/test
test_mi_graph <- qgraph(test_mi, layout='spring', vsize=2)

##heatmap

heatmap(as.matrix(test))

#better from gplots
my_col_breaks <- c(0,1,2,3,4)

my_pallete <- colorRampPalette(c("red", "white"))((n=4))

heatmap.2(as.matrix(test), scale = "none", trace = "none", symbreaks = F, col = my_pallete, symm = F, 
          breaks = my_col_breaks, margins=c(6,6), cexRow = .1, cexCol = .1,
          hclustfun = function(x) hclust(x, method="ward.D2")
          )

?hclust

?heatmap.2

## new analysis with global alignment
test <- stringDist(u_v4, method = "levenshtein", ignoreCase = TRUE, 
                   diag = FALSE,
                   upper = FALSE,
                   type = "global")
View(as.data.frame(test))

u_v4[grep("IGHV5-15", names(u_v4), value = T)]

test <- stringDist(u_v4[grep("IGHV5-15", names(u_v4), value = T)], method = "levenshtein", ignoreCase = TRUE, 
          diag = FALSE,
          upper = FALSE,
          type = "overlap",
          )

heatmap.2(as.matrix(test), scale = "none", trace = "none", symbreaks = F, symm = F, 
                  margins=c(6,6))

plot(h1)
"IGHV5-157*02_S3096"
?pairwiseAlignment

IGHV4-NGC42_S4754
IGHV4-72*01 IGHV4-NGC35 IGHV4-67*01_S0428 IGHV4-149*02_S9930

u_v4["IGHV5-157*02_S3096"]
pair_wise <- pairwiseAlignment(pattern = u_v4["IGHV4-NGC38_S3447"],
                  subject = u_v4["IGHV4-NGC38"],
                  type = "global")
pair_wise
p2 <- alignedSubject(pair_wise)
p3 <- alignedPattern(pair_wise)


g2 <- ggmsa(c(p2,p3), color = 'Zappo_AA',
            seq_name = T, font = 2)

print(g2)

library(DECIPHER)
seq <- c(aligned(pattern(pair_wise)), aligned(subject(pair_wise)))
BrowseSeqs(seq)

##DIST CONVERSION FUNCTION TO DATAFRAME
.myfun <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B)-1):1),
    value = as.vector(inDist))
}
x <- .myfun(test)
grep


## msa

sel <- c("IGHV3-53*01_S2000",
         "IGHV3-53*01",
         "IGHV3-53*01_S4256",
         "IGHV3-53*01_S3647",
         "IGHV3-122*01",
         "IGHV3-122*01_S8388",
         "IGHV3-122*01_S8946",
         "IGHV3-122*01_S0689",
         "IGHV3-53*02_S5299",
         "IGHV3-53*02")
sel <- u_v4[sel]

library(msa)
sel_alg <- msaMuscle(sel, type = "dna",
                     verbose = TRUE)
sel_alg_seq <- DNAStringSet(sel_alg)
DECIPHER::BrowseSeqs(sel_alg_seq)
writeXStringSet(sel, "~/Box Sync/rsv-sc/results/2021-02-10/sequences-test.fasta")

##plot only the ones with less than 4 similarities

x_dim <- x %>% filter(value <= 3) 
x_dim_dcast <- reshape2::dcast(x_dim, row ~ col)
rownames(x_dim_dcast) <- x_dim_dcast$row
x_dim_dcast <- as.matrix(x_dim_dcast[-1])
is.na(x_dim_dcast)

##renaming
x_dim_dcast[is.na(x_dim_dcast)] <- 0

heatmap.2(x_dim_dcast, scale = "none", trace = "none", symbreaks = F, symm = F, 
          margins=c(6,6))
