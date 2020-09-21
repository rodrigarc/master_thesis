BiocManager::install("PCAtools")
install.packages("factoextra")
install.packages("Seurat")
install.packages("igraph")
install.packages("png")
install.packages("hdf5r")
library(hdf5r)
library(png)
library(igraph)
library(Seurat)
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
library(M3C)



loadfonts(device = "win")
devtools::install_github("r-lib/devtools")
devtools::install_github("richardjtelford/ggbiplot", ref = "experimental")

#data processing, combining dataframes etc.
sf_10 <- sf_trim_0.1$summaries %>% 
  mutate(sf = 10, 
         trimmed.integer = as.integer(trimmed.mean.quality)) 
sf_10[is.na(sf_10)] <- 0


sf_13 <- sf_trim_0.0501187233627272$summaries %>% 
  mutate(sf = 13, 
         trimmed.integer = as.integer(trimmed.mean.quality))
sf_13[is.na(sf_13)] <- 0


sf_16 <- sf_trim_0.0251188643150958$summaries %>% 
  mutate(sf = 16, 
         trimmed.integer = as.integer(trimmed.mean.quality))
sf_16[is.na(sf_16)] <- 0


sf_20 <- sf_trim_0.01$summaries %>% 
  mutate(sf = 20, 
         trimmed.integer = as.integer(trimmed.mean.quality)) 
sf_20[is.na(sf_20)] <- 0


sf30 <- sf_30$summaries %>% 
  mutate(sf = 30, 
         trimmed.integer = as.integer(trimmed.mean.quality)) 
sf30[is.na(sf30)] <- 0

sf40 <- sf_40$summaries %>% 
  mutate(sf = 40, 
         trimmed.integer = as.integer(trimmed.mean.quality)) 
sf40[is.na(sf40)] <- 0

sf50 <- sf_50$summaries %>% 
  mutate(sf = 50, 
         trimmed.integer = as.integer(trimmed.mean.quality))
sf50[is.na(sf50)] <- 0

combined <- rbind(sf_10,sf_13,sf_16,sf_20,sf30,sf40,sf50)
combined <- combined[,c(2:14)]

#rename for better PCA anmes

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

#do pca based on filtered vs non-filtered
test_f <- test %>% mutate(filtering = ifelse(
  test$file_ID %in% sf_filtered$file_ID, "Filtered", "Not filtered"))


pca <- PCA(combined_umap[,c(-1,-2,-3,-13,-14)])
fviz_pca_var(pca, labelsize =5, repel = T, axes.linetype = "dotted", col.circle = "white")+
  ggtitle("")+
  theme_cowplot()+
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 16, face = "bold", color = "black"),
        panel.grid = element_line("white"),
        axis.ticks = element_line(size = 1.5),
        axis.line = element_line(color = "black", size = 1.5))


fviz_pca_ind(pca, label = "var", 
                habillage = as.factor(combined_umap$filtered),
                palette = c("#440154ff","#DCE319FF"),
                alpha.ind = 0.7,
             title = "", addEllipses = F)+
  theme_cowplot() +
  theme(legend.position = "right", 
        legend.direction = "vertical",
        legend.title = element_blank(), 
        axis.line = element_line(color = "black", size = 1.5),
        axis.ticks = element_line(size = 1.5),
        axis.text = element_text(face = "bold", size = 18),
        axis.title = element_text(face ="bold"),
        rect = element_rect(fill = "transparent"),
        text = element_text(size = 20))
                

#try UMAP
combined_umap <- sf_20 %>% mutate(
  filtered = ifelse(trimmed.mean.quality > 30, "> 30", "< 30")) %>%
  select(-c(folder.name, file.name, sf, file.path,trimmed.integer, filtered))

combined_umap <- sf_20 %>% 
  filter(trimmed.mean.quality >30 & raw.length > 400) %>% 
  select(-c(folder.name, file.name, sf, file.path,trimmed.integer)) 
combined_umap <- combined_umap %>% select(raw.length, trim.start, trim.finish, 
                                          trimmed.mean.quality,
                                          trimmed.secondary.peaks)

kmeans(combined_umap[,c(-6)], 3)

if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
    Sys.info()["sysname"] == "Darwin" && getRversion() == "4.0.0") {
  parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
}

#testing methods to calculate the best klustering number
combined_m3c <- M3C((combined_umap[,c(-6)]))

elbow <- fviz_nbclust(combined_umap[,c(-6)], kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method")
elbow
silhouette <- fviz_nbclust(combined_umap[,c(-6)], kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
silhouette
gap_stat <- fviz_nbclust(combined_umap[,c(-6)], kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")
gap_stat

M3C::umap(t(combined_umap[,-11]),
     labels=as.factor(t(combined_umap$cluster)), 
     dotsize = 0.5)

#test another pca

pca <- PCAtools::pca(combined_umap[,c(-1,-2,-3,-13,-14)], removeVar = 0.1)
screeplot(pca)
biplot(pca)
pairsplot(pca)
horn <- parallelPCA()
elbow <- findElbowPoint(pca$variance)
screeplot(pca,
          components = getComponents(pca, 1:9),
          vline = c(elbow)) +
  geom_text(aes(elbow + 1, 50, label = "Elbow", vjust = -1))

biplot(pca)

#test unsupervised clustering
test <- kmeans(m, 3)
combined_umap$cluster <- test$cluster
test

pca_x <- PCAtools::pca(combined_umap[,-11])
PCAtools::biplot(pca_x)
hist(combined_umap,
    breaks = 200)

library("ggpubr")
ggdensity(combined_umap$raw.length, fill = "lightgray") #need to normalize
ggqqplot(combined_umap$trimmed.mean.quality)

#test normalization
m <- sweep(combined_umap, 2, colSums(combined_umap), FUN="/")
m <- scale(m, center=T, scale=colSums(m))
library(ggpubr)
ggdensity(m[,4], fill = "lightgray") #need to normalize
ggqqplot(m[,4])

pca <- PCAtools::pca(m)
screeplot(m_pca)
library(factoextra)
library(PCAtools)
pca_summary <- FactoMineR::summary.PCA(pca)
pca <- pca(m)
screeplot(pca)
biplot(pca)
horn <- parallelPCA(m)
horn$n
elbow <- findElbowPoint(pca$variance)
elbow
library(ggplot2)

screeplot(pca,
          components = getComponents(pca, 1:8),
          vline = c(horn$n, elbow)) +
  geom_text(aes(horn$n + 1, 50, label = "Horn's", vjust = -1)) +
  geom_text(aes(elbow + 1, 50, label = "Elbow", vjust = -1))
plotloadings(pca,
             rangeRetain = 0.01,
             labSize = 3.0,
             title = 'Loadings plot',
             subtitle = 'PC1, PC2, PC3, PC4, PC5',
             caption = 'Top 1% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

library("FactoMineR")
res.pca <- PCA(m, graph = FALSE)
var <- get_pca_var(res.pca)
library("corrplot")
corrplot(var$contrib, is.corr=FALSE)  
fviz_cos2(res.pca, choice = "var", axes = 1:2)
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)




M3C::umap(t(m),
          dotsize = 0.5)

gap_stat <- factoextra::fviz_nbclust(m, kmeans, nstart = 25,  k.max = 15,
                                     method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")
gap_stat

#testing hierarchical clustering
##cluster normalized - hierarchical cluster
summary(m)
dist_mat <- dist(m, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 30)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 30)
library(dendextend)
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, h = 30)
plot(avg_col_dend)
seeds_df_cl <- mutate(data.frame(m), cluster = cut_avg)
g5 <- ggplot(seeds_df_cl, aes(x=trimmed.mean.quality, 
                        y = raw.length, color = factor(cluster))) + 
  geom_point()
g5
g6 <- umap(t(seeds_df_cl[,-11]),
labels=as.factor(seeds_df_cl$cluster),
dotsize = 0.5)
g6
grid.arrange(g5, g6, nrow = 1)

library(kableExtra)
seeds_df_cl %>% group_by(cluster) %>% summarise_all(mean) %>% kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

#cluster normalized - kmeans
m <- sweep(combined_umap, 2, colSums(combined_umap), FUN="/")
m <- scale(m, center=T, scale=colSums(m))
cluster<- kmeans(m,3)
m_df <- as.data.frame(m)
m_df$cluster <- cluster$cluster

g3<- ggplot(m_df, aes(x=trimmed.mean.quality, 
                        y = raw.length, color = factor(cluster))) + 
  geom_point()
g3
g4 <- umap(t(m_df[,-11]),
               labels=as.factor(m_df$cluster),
               dotsize = 0.5)
g4
library(gridExtra)
grid.arrange(g3, g4, nrow = 1)

library(kableExtra)
m_df %>% group_by(cluster) %>% summarise_all(mean) %>% kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F)

#################cluster not normalized - kmeans
cluster<- kmeans(combined_umap,3)
combined_umap$cluster <- cluster$cluster

library(M3C)
library(ggplot2)
g2 <- umap(t(combined_umap[,-11]),
     labels=as.factor(combined_umap$cluster),
     dotsize = 0.5)
g2
g1 <- ggplot(combined_umap, aes(x=trimmed.mean.quality, 
                 y = raw.length, color = factor(cluster))) + 
  geom_point()
g1


library(gridExtra)
grid.arrange(g1, g2, nrow = 1)

library(kableExtra)
library(xlsx)
library(tidyverse)
combined_umap %>% group_by(cluster) %>% summarise_all(sd) %>% kable() %>%
  kable_styling(bootstrap_options = "striped", full_width = F)


combined_umap %>% group_by(cluster) %>% summarise_all(funs(mean,sd, n())) %>% 
  write.csv("~/Desktop/rodrigo/R filtering improvement/2020-07-15_cluster-kmeans-not-normalized_mean-sd-count.csv")

apply(asplit(my_array, 1:2), 2, unlist)
sf_mean <- mean(sf_filtered$raw.length)
hist(sf_filtered$raw.length,breaks = 100, prob = T)

m<-mean(sf_filtered$trim.finish)
std<-sqrt(var(sf_filtered$trim.finish))
hist(sf_filtered$trim.finish, prob=T,main="Trim start", breaks = 100)
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE)

hist(combined_umap$trim.start, breaks = 100)

res <- M3C(combined_umap, removeplots = TRUE, iters=2,
           objective='PAC', fsize=8, lthick=1, dotsize=1.25, cores = 1)

#testing another algorith for clustering
install.packages("dbscan")
library(dbscan)

cl <- dbscan(combined_umap[,-11], eps = 2 )
plot(combined_umap, col=cl$cluster+1, pch=20)
combined_umap$cluster <- cl$cluster
ggplot(combined_umap, aes(x=trimmed.mean.quality, 
                          y = raw.length, color = factor(cluster))) + 
  geom_point()
umap(t(combined_umap[,-11]),
     labels=as.factor(combined_umap$cluster),
     dotsize = 0.5)
