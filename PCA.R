BiocManager::install("PCAtools")
install.packages("factoextra")
library(openxlsx)
library(reshape2)
library(plyr)
library(scales)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(devtools)
library(dplyr)
library(grid)
library(pkgbuild)
library(devtools)
library(PCAtools)
library(FactoMineR)
library(factoextra)
library(extrafont)

loadfonts(device = "win")
devtools::install_github("r-lib/devtools")
devtools::install_github("richardjtelford/ggbiplot", ref = "experimental")

test <- sf$summaries %>% select(trimmed.mean.quality, 
                                raw.mean.quality,
                                trimmed.secondary.peaks,
                                trim.finish,
                                trim.start,
                                raw.length,file_ID) %>%
  rename( "trimmed.mean.quality" = "T Mean quality", 
         "raw.length" = "Length" ,
         "trim.finish"="Trim end" ,
         "trim.start" = "Trim start" ,
         "trimmed.secondary.peaks" = "T sec. peaks",
         "raw.mean.quality" = "Mean quality")

test_f <- test %>% mutate(filtering = ifelse(
  test$file_ID %in% sf_filtered$file_ID, "Filtered", "Not filtered"))
pca <- PCA(test_f[,c(-7,-8)])
fviz_pca_var(pca, labelsize =5, repel = T, axes.linetype = "dotted", col.circle = "white")+
  xlab("PC1 (58%)") + 
  ylab("PC2 (20.5%)")+ 
  ggtitle("")+
  theme_cowplot()+
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, face = "bold", color = "black"),
        panel.grid = element_line("white"),
        axis.ticks = element_line(size = 1.5),
        axis.line = element_line(color = "black", size = 1.5))

fviz_pca_ind(pca, label = "var", 
                habillage = as.factor(test_f$filtering),
                palette = c("#440154ff","#DCE319FF"),
                alpha.ind = 0.7,
             title = "", addEllipses = F)+
  theme_cowplot() + xlab ("PC1 (58%)") + ylab("PC2 (20.5%)")+
  theme(legend.position = "right", 
        legend.direction = "vertical",
        legend.title = element_blank(), 
        axis.line = element_line(color = "black", size = 1.5),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face ="bold"),
        rect = element_rect(fill = "transparent"),
        text = element_text(size = 20))
                