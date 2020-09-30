#Sys.setenv(LANG = "en")
#Sys.setlocale("LC_ALL", "C")
library(Biostrings)
library(msa)
library(ape)
library(ggtree)
library(tidyverse)
library(ggrepel)
library(gridExtra)
setwd("~/Box Sync/P6-Germlines/Databases - FASTA files/Vasquez-Bernat")
wd <- getwd()

#load fasta file
bernat_db <- readDNAStringSet("1198.fasta",
                             format="fasta")
    names(bernat_db) <- sub("IG","vb_IG", names(bernat_db))

    bernat_db <- unique(bernat_db)
NGCdb <- readDNAStringSet("~/Box Sync/Rodrigo & Klara & Sebastian/Germlines/NGC candidates/NGCs_v7.fasta",
                                         format="fasta")
    names(NGCdb) <- sub("IG","ngc_IG", names(NGCdb))
          
#Combine all individualized databases
combined <- c(NGCdb,bernat_db)

aligned_Muscle <- msaMuscle(combined, verbose = TRUE)
#aligned_Omega <- msaClustalOmega(combined, verbose = TRUE)
#aligned_W <- msaClustalW(combined, verbose = TRUE)

aligned_bin <- as.DNAbin(aligned_Muscle)
d <- dist.dna(aligned_bin)

tree <- nj(d)
g <- ggtree(tree, layout = "rectangular", MAX_COUNT = 1)
g+geom_tiplab(size = 0.8)
g1 <- g$data %>% mutate(label2 = ifelse(grepl("vb_", label),"bernat",
                                 ifelse(grepl("ngc_", label), "lore",
                                 "")))

g2 <- ggtree(g1, aes(color = label2)) +
  geom_tippoint(size = 0) + 
  ggtitle("Bernat db vs Loré db")+
  theme(legend.position="top",legend.text=element_text(size=25), 
        legend.key.size = unit(6,"line"),
        plot.title = element_text(size = 30, face = "bold")) + 
  labs(color = "")+
  scale_color_manual(values = c("black","red","blue","black"))+
  geom_tiplab(size = .5, hjust = 0.1)+
  geom_treescale(width = 0.1, x = .4)
g2

#ggsave("db_v4_Muscle.pdf",
       width = 40, height = 40, 
       units = "cm", limitsize = FALSE)

