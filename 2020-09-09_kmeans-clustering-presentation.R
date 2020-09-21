#####################################################################################
################## K-MEANS CLUSTERING ###############################################
library(M3C)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(factoextra)
library(stringr)
library(Biostrings)

#read sf_trim 10-20 file to be faster (~/Box Sync/codes/trims)
#UMAP
load("~/Box Sync/codes/trims/sf_trim_10_20.RData")
sf_20 <- sf_trim_0.01$summaries %>%
  mutate(plate = paste(substring(folder.name,5,6)))%>%
  mutate(well = gsub("\\-.*|\\_.*", "", file.name)) %>%
  mutate(ID = substring(folder.name, 1, 3)) %>%
  mutate(chain = "HC")%>% 
  mutate(file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))%>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()
sf_20[is.na(sf_20)] <- 0

#read df with cdr3
df <- read.csv("~/Box Sync/codes/R filtering improvement/CDR3_sec_pek_df.csv")
sf_20_cdr3<- merge(sf_20,df)

#edditing file without some columns and filtering > 400 length
combined_umap <- sf_20_cdr3 %>% 
  filter(raw.length > 400) %>% 
  select(-c(folder.name, file.name,file.path,
             well, plate, ID, chain)) %>%
  relocate(file_ID)
  

##testing the best number of clusters
elbow <- fviz_nbclust(combined_umap, kmeans, method = "wss") +
  labs(subtitle = "Elbow method")
elbow + geom_vline(xintercept = 3)


#plotting umap colorring by filters
#use the R filtering data improvement to combine the CDR3 data fasta
sf_20_filtering <- sf_20_cdr3 %>% mutate(f4  = ifelse(raw.length >= 400 &
                                                  trim.start < 33 &
                                                  trim.finish > 409 &
                                                  trimmed.mean.quality >= 30,
                                                1 , 0),
                                  f4_peak = ifelse(raw.length >= 400 &
                                                       trim.start < 33 &
                                                       trim.finish > 409 &
                                                       trimmed.mean.quality >= 30 &
                                                       trimmed.secondary.peaks <= 50,
                                                     1 , 0),
                                  f4_peak_cdr3 = ifelse(raw.length >= 400 &
                                                            trim.start < 33 &
                                                            trim.finish > 409 &
                                                            trimmed.mean.quality >= 30 &
                                                            trimmed.secondary.peaks <= 50 &
                                                            sec.peak.CDR3 <= 5,
                                                          1 , 0))

sf_20_filters <- sf_20_filtering %>%
  mutate(Total = select(., c(f4, f4_peak, f4_peak_cdr3)) %>% rowSums(na.rm = TRUE)) %>%
  select(-c(f4,f4_peak, f4_peak_cdr3, file.name, folder.name, file.path,  
             well, ID, chain, plate, file_ID))
#replacing 0,1,2,3 to Filtered out, 4f , 4f_peak, 4f_peak_cdr3
sf_20_filters$Total[sf_20_filters$Total == 0] <- "Filtered out = 4F"
sf_20_filters$Total[sf_20_filters$Total == 1] <- "Filtered out = 4F + P "
sf_20_filters$Total[sf_20_filters$Total == 2] <- "Filtered out = 4F + P + C"
sf_20_filters$Total[sf_20_filters$Total == 3] <- "Filtered in"
sf_20_filters$sec.peak.CDR3 <- as.numeric(sf_20_filters$sec.peak.CDR3)

g4 <- umap(t(sf_20_filters[,-12]),
           labels=as.factor(sf_20_filters$Total),
           dotsize = 0.5, seed = 250, 
           scale = 3,
           colvec = c("#33CCFF", "black", "yellow", "red"),
           controlscale = T)
g4 <- g4 + labs(x = "UMAP1", y = "UMAP2")
g4


#clustering
set.seed(250)
cluster<- kmeans(sf_20_filters[,-12], 4)
sf_20_filters$cluster <- cluster$cluster

#plotting umap colorring by kmeans
g2 <- umap(t(sf_20_filters[,-c(12,13)]),
           labels=as.factor(sf_20_filters$cluster),
           dotsize = 0.5, seed = 250, printwidth = 20
)
g2 <- g2 + labs(x = "UMAP1", y = "UMAP2")
g2

g1 <- ggplot(combined_umap, aes(x=trimmed.mean.quality, 
                                y = raw.length, color = factor(cluster))) + 
  geom_point()

g1
grid.arrange(g1, g2, nrow = 1)

###plotting UMAP old filter
sf_20_old_dedup<- sf_20_cdr3 %>% mutate(filtered  = ifelse(raw.length > 300 & 
                                                          raw.length < 600 &
                                                          trim.start < 65 & 
                                                          trim.finish > 150 &
                                                          trimmed.mean.quality > 30 & 
                                                          trimmed.secondary.peaks < 150,
                                                          "Filtered in" , "Filtered out"))
sf_20_old_dedup_umap <- sf_20_old_dedup %>%
  select(-c(file.name, folder.name, file.path,  
            well, ID, chain, plate, file_ID))
                                         
g5 <- umap(t(sf_20_old_dedup_umap[,-12]),
           labels=as.factor(sf_20_old_dedup_umap$filtered),
           dotsize = 0.5, seed = 250,
           printwidth = 20)
g5 <- g5 + labs(x = "UMAP1", y = "UMAP2")

###plotting umap 4fv2 - changing <32 to <50 trim.start
sf_20_filtering_4f_v2 <- sf_20_cdr3 %>% mutate(f4  = ifelse(raw.length >= 400 &
                                                        trim.start < 50 &
                                                        trim.finish > 409 &
                                                        trimmed.mean.quality >= 30,
                                                      1 , 0),
                                         f4_peak = ifelse(raw.length >= 400 &
                                                            trim.start < 50 &
                                                            trim.finish > 409 &
                                                            trimmed.mean.quality >= 30 &
                                                            trimmed.secondary.peaks <= 50,
                                                          1 , 0),
                                         f4_peak_cdr3 = ifelse(raw.length >= 400 &
                                                                 trim.start < 50 &
                                                                 trim.finish > 409 &
                                                                 trimmed.mean.quality >= 30 &
                                                                 trimmed.secondary.peaks <= 50 &
                                                                 sec.peak.CDR3 <= 5,
                                                               1 , 0))

sf_20_filters_4f_v2 <- sf_20_filtering_4f_v2 %>%
  mutate(Total = select(., c(f4, f4_peak, f4_peak_cdr3)) %>% rowSums(na.rm = TRUE)) %>%
  select(-c(f4,f4_peak, f4_peak_cdr3, file.name, folder.name, file.path,  
            well, ID, chain, plate, file_ID))
#replacing 0,1,2,3 to Filtered out, 4f , 4f_peak, 4f_peak_cdr3
sf_20_filters_4f_v2$Total[sf_20_filters_4f_v2$Total == 0] <- "Filtered out = 4F"
sf_20_filters_4f_v2$Total[sf_20_filters_4f_v2$Total == 1] <- "Filtered out = 4F + P "
sf_20_filters_4f_v2$Total[sf_20_filters_4f_v2$Total == 2] <- "Filtered out = 4F + P + C"
sf_20_filters_4f_v2$Total[sf_20_filters_4f_v2$Total == 3] <- "Filtered in"
sf_20_filters_4f_v2$sec.peak.CDR3 <- as.numeric(sf_20_filters_4f_v2$sec.peak.CDR3)

g6 <- umap(t(sf_20_filters_4f_v2[,-12]),
           labels=as.factor(sf_20_filters_4f_v2$Total),
           dotsize = 0.5, seed = 250, 
           scale = 3,
           colvec = c("#33CCFF", "black", "yellow", "red"),
           controlscale = T,
           printwidth = 20)
g6 <- g6 + labs(x = "UMAP1", y = "UMAP2")
g6
##################################################################################
#writing CSV containing the mean, sd and coutns of each cluster
x <- sf_20_filters %>% group_by(cluster) %>% summarise_all(funs(mean,sd, n()))
write.csv(x ,
            "~/Box Sync/Gunilla presentation/2020-09-09_cluster-kmeans-not-normalized_mean-sd-count_2.csv")

#saving UMAP as TIFF 
##save UMAP filters 4f
ggsave(g4, filename = "~/Box Sync/Gunilla presentation/UMAP_filters.tiff",
       limitsize = FALSE, bg = "transparent", device = "tiff", width = 23,
       height = 14, dpi = 150, units = "cm")

#save old_filter
ggsave(g5, filename = "~/Box Sync/Gunilla presentation/UMAP_old_filter.tiff",
       limitsize = FALSE, bg = "transparent", device = "tiff", width = 14,
       height = 9, dpi = 150, units = "cm")

#save kmeans
ggsave(g5, filename = "~/Box Sync/Gunilla presentation/UMAP_Kmeans.tiff",
       limitsize = FALSE, bg = "transparent", device = "tiff", width = 12,
       height = 9, dpi = 150, units = "cm")

