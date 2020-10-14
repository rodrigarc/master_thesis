library(ggtree)
library(treeio)
library(ggplot2)
library(Biostrings)

###### making fasta names unique ######
setwd("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results")
ls_fa <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-14/msa",
                    full.names = T)
ls_fa <- grep(".fasta", ls_fa, value = T)
ls_fa_names <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-14/msa")
ls_fa_names <- grep(".fasta", ls_fa_names, value = T)
fa <- lapply(ls_fa, readBStringSet)
names(fa) <- ls_fa_names
for (i in names(fa)){
  names(fa[[i]]) <- make.unique(names(fa[[i]]))
  writeXStringSet(fa[[i]],paste0("2020-10-14/msa",i))
}


##### doing fasttree using the bash code in src #####


####plotting trees ####
ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-14/trees/", 
                 full.names = T)
names <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-14/trees/", 
                    full.names = F)
names <- grep(".tre",names,value = T, fixed = T)

ls <- grep(".tre",ls,value = T, fixed = T)
trees <- lapply(ls, read.tree)
names(trees) <- names

trees_gg <- lapply(trees, function(x) ggtree (x) + geom_tiplab(size = 0.2))

for (i in names(trees_gg)){
  trees_gg[[i]] + ggtitle(grep(i, names(trees_gg), value = T)) + 
  ggsave(paste0("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-14/plot/",i,".pdf"),
         width = 40, height = 40, units = "cm")
}
