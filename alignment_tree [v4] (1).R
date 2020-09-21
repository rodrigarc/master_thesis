#Sys.setenv(LANG = "en")
#Sys.setlocale("LC_ALL", "C")
library(Biostrings)
library(msa)
library(ape)
library(ggtree)
library(tidyverse)
library(ggrepel)
library(gridExtra)
setwd("~/Box Sync/Students/Rodrigo & Klara & Sebastian/Germlines/")
wd <- getwd()

#load fasta file
watsondb <- readDNAStringSet("~/Box Sync/Students/Rodrigo & Klara & Sebastian/Germlines/watson/V.fasta",
                             format="fasta")
            names(watsondb) <- sub("IG","ws_IG", names(watsondb))

NGCdb <- readDNAStringSet("~/Box Sync/Students/Rodrigo & Klara & Sebastian/Germlines/NGC candidates/NGCs_v7.fasta",
                                         format="fasta")
            names(NGCdb) <- sub("IG","ngc_IG", names(NGCdb))
            
db_v4 <- readDNAStringSet("~/Box Sync/Students/Rodrigo & Klara & Sebastian/Germlines/IgM_v4/germlines_v4_unique.fasta",
                                format="fasta")
            names(db_v4) <- sub("IG","v4_IG", names(db_v4))

#Combine all individualized databases
combined <- c(db_v4,NGCdb,watsondb)

aligned_Muscle <- msaMuscle(combined, verbose = TRUE)
aligned_Omega <- msaClustalOmega(combined, verbose = TRUE)
#aligned_W <- msaClustalW(combined, verbose = TRUE)

aligned_bin <- as.DNAbin(aligned_Muscle)
d <- dist.dna(aligned_bin)

tree <- nj(d)
g <- ggtree(tree, layout = "rectangular", MAX_COUNT = 1)
g+geom_tiplab(size = 0.8)
g1 <- g$data %>% mutate(label2 = ifelse(grepl("ws_", label),"Watson",
                                 ifelse(grepl("ngc_", label), "NGC",
                                 ifelse(grepl("v4_", label), "ID_v4",""))))

g2 <- ggtree(g1, aes(color = label2)) +
  geom_tippoint(size = .1) + 
  ggtitle("Individualized v4 vs Watson/NGC")+
  theme(legend.position="top",legend.text=element_text(size=25), 
        legend.key.size = unit(6,"line"),
        plot.title = element_text(size = 30, face = "bold")) + 
  labs(color = "")+
  scale_color_manual(values = c("black","red","blue","black"))+
  geom_tiplab(size = .8, hjust = 0.1)+
  geom_treescale(width = 0.1, x = .4)
g2
ggsave("db_v4_Muscle.pdf",
       width = 40, height = 40, 
       units = "cm", limitsize = FALSE)

