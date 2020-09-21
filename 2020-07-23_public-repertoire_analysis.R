install.packages("ggnewscale")

library(data.table)
library(kableExtra)
library(msa)
library(Biostrings)
library(dplyr)
library(ape)
library(ggtree)
library(ggplot2)
library(ggnewscale)

setwd("~/Desktop/rodrigo/public-clonotype/Analysis/")
#selecting files of clonoquery and reading it
ls <- list.files("~/Desktop/rodrigo/public-clonotype/Public_clonoquery", full.names = T)
ls <- grep("full", ls, value = T)
tab <- lapply(ls, fread, fill = T, blank.lines.skip = T, 
               stringsAsFactors = F,
               header = T)
#remove empty rows with 0 clones 
tab_n <- lapply(tab, na.omit)

#combine list into df
tab_c <- rbindlist(tab_n, fill = T)

#test some group_by
tab_test <-  tab_c %>% group_by(V_gene, J_gene, CDR3_length) %>% 
  summarise(ID = paste(unique(name), collapse=", "), counts = sum(count))
kable(tab_test) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


#alignment and tree
aa <- Biostrings::AAStringSet(tab_c$CDR3_aa)
aa_vdj <- Biostrings::AAStringSet(tab_c$VDJ_aa)

  
#trying aligning with ugpma and nj, nj it is better for evolutionary studies
aligned_Muscle <- msaMuscle(aa, verbose = T)
aligned_vdj_Muscle <- msaMuscle(aa_vdj, verbose = TRUE)
aligned_muscle_nj <- msaMuscle(aa, verbose = TRUE, cluster = "neighborjoining")
aligned_vdj_Muscle


#calculate distance and plot tree
aligned_bin <- as.AAbin(aligned_Muscle)
d <- dist.aa(aligned_bin)
tree <- nj(d)
g <- ggtree(tree, layout = "circular",branch.length = "none", MAX_COUNT = 1)
g


#write multiple sequence alignment
write.FASTA(aligned_bin, 
          "~/Desktop/rodrigo/public-clonotype/Analysis/2020-07-24_aa-vdj-aligned.fasta")


#remember to annotate the tree
public <- c("IGHV4-117","IGHV3-50",
            "IGHV4-NGC30","IGHV4-NGC23", 
            "IGHV3-NGC11","IGHV4-NGC42","IGHV4-NGC44")

tab_cl <- tab_c %>% mutate(node = 1:nrow(tab_c), 
                           ID = name,
                           "V gene family" = ifelse(V_gene %in% public,
                                                    public,substring(V_gene,1,5))) %>% 
  select(node, ID, everything())

#only IGHV4-NGC30 & IGHJ4-3
tab_cl <- tab_c %>% mutate(node = 1:nrow(tab_c), 
                           ID = name,
                           "V gene family" = ifelse(V_gene == "IGHV4-NGC30" & J_gene == "IGHJ4-3",
                                                    "IGHV4-NGC30 x IGHJ4-3","")) %>% 
  select(node, ID, everything())

p <- g %<+% tab_cl 
p1 <- p + geom_tippoint(aes(color = ID)) + 
    theme(legend.position = "right")


p2 <- p + new_scale_fill()
p3 <- gheatmap(p2, tab_cl[,c(2)], offset=.8, width=.1,
               colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_brewer(palette = "Set1", name = "ID")

p3 <- p3 + new_scale_fill()
p4 <- gheatmap(p3, tab_cl[,c(8)], offset=15, width=.2,
               colnames_angle=90, colnames_offset_y = .25) + 
  scale_fill_viridis_c(name = "CDR3 Length")

p4 <- p4 + new_scale_fill()
p5 <- gheatmap(p3, tab_cl[,c(30)], offset=30, width=.3,
               colnames_angle=90, colnames_offset_y = .25) + 
  scale_fill_brewer(palette = "Set3", name = "Public clonotype\n(> 3 animals)")
p5
p5 <- p5 + new_scale_fill()
p5
p6 <- gheatmap(p5, tab_cl[,c(7)], offset=200, width=.3,
                     colnames_angle=90, colnames_offset_y = .25) + 
  scale_fill_viridis_d(option = "plasma", name = "J genes")
p6

nrow()levels()msaplot(g, "2020-07-24_aa-cdr3-aligned.fasta", bg_line = F)
  

  
  
  
  