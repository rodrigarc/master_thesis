install.packages("pheatmap")
library(pheatmap)

protein_sequences <- system.file("extdata", "sample.fasta", package = "ggmsa")
aln <- readAAMultipleAlignment("_out.200528181331420pz1eufv2MBodUTkthqZLGlsfnormal.fasta")
ggmsa(protein_sequences, start = 1, end = 300, color = "Chemistry_AA",font = NULL)
protein_sequences
aln
aln <- readAAMultipleAlignment("~/Desktop/vinnycius/msa_renamed.fasta")
aln <- unmasked(aln)
names(aln)[1]
ref <- aln[1]
bm <- sapply(1:length(aln),function(i){
  as.numeric(as.matrix(aln[i])==as.matrix(ref))
})
bm <- t(bm)
rownames(bm) <- names(aln)
heatmap <- pheatmap(bm[nrow(bm):1,1:1273],cluster_rows=T,cluster_cols=F, 
                    color = c("#FF6666","#CCE5FF"), legend_breaks = c(0,1),
                    legend_labels = c("Divergent", "Convergent"))
heatmap(bm)         
