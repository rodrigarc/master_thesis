install.packages("pheatmap")
library(pheatmap)

setwd("/Users/karlor/Desktop/sebastian/")

aln <- readAAMultipleAlignment("msa_renamed_final.fasta")
aln <- unmasked(aln)
names(aln)[1]
ref <- aln[1]
bm <- sapply(1:length(aln),function(i){
  as.numeric(as.matrix(aln[i])==as.matrix(ref))
})
bm <- t(bm)
rownames(bm) <- names(aln)
heatmap <- pheatmap(bm[nrow(bm):1,1:1273],cluster_rows=T,cluster_cols=F, 
                    color = c("#FF6666","#F8F8F8"), legend_breaks = c(0:1), 
                    cellwidth = 0.4, cellheight = 15, annotation_names_col = T,
                    legend_labels = c("Divergent", "Convergent"))
       
