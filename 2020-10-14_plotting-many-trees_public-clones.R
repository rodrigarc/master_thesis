library(ggtree)
library(treeio)
library(ggplot2)
library(Biostrings)
library(dplyr)


###### making fasta names unique ######
###### at the end I did not use that, I just created unique names before doing the alignment ######
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
ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-15/trees/", 
                 full.names = T)
names <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-15/trees/", 
                    full.names = F)
names <- grep(".tre",names,value = T, fixed = T)

ls <- grep(".tre",ls,value = T, fixed = T)
trees <- lapply(ls, read.tree)
names(trees) <- names

trees_gg <- lapply(trees, ggtree)
trees[[1]]
trees_gg[[3]] + geom_b
for(i in names(trees_gg)){
  trees_gg[[i]] <- trees_gg[[i]]$data %>% mutate(label2 = ifelse(grepl("B2_", label),"B2",
                                                                 ifelse(grepl("B1_", label), "B1",
                                                                        ifelse(grepl("PV_", label), "PV",""))),
                                                 label3 = substring(label, 7,9))

  }
root(trees[[1]], 
     node = trees[[1]][["parent"]][trees[[1]][["label3"]] == "con"], 
     outgroup = "testing")

trees_gg[[1]][["label"]][trees_gg[[1]][["label3"]] == "con"]

for(i in names(trees_gg)){
  trees_gg[[i]] <- as.treedata(trees_gg[[i]])
  trees_gg[[i]] <- root(trees_gg[[i]], 
                        node = trees_gg[[i]][["parent"]][trees_gg[[i]][["label3"]] == "con"])
}
 
trees_gg <- lapply(trees_gg, ggtree)

for (i in names(trees_gg)){
  trees_gg[[i]] + geom_tiplab(size = 0.2) + 
    geom_tippoint(aes(color=label3, shape = label2), size = 2)+
    ggtitle(grep(i, names(trees_gg), value = T)) + 
    geom_treescale()
  ggsave(paste0("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-15/plot/",i,".pdf"),
         width = 40, height = 40, units = "cm")
}

