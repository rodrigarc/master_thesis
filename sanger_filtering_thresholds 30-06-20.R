library(dplyr)
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(tidyverse)


#sanger_filtering_thresholds_tests

sf <- summarise.abi.folder("~/Desktop/rodrigo/sanger_HC/sanger_test",                           secondary.peak.ratio = 0.33, 
                           trim.cutoff = 0.01)
View(sf)

#create same columns for the data frame not filtered
sf_summary <- sf_40$summaries %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

sf_summary[is.na(sf_summary)] <- 0

#test only here
sf_filtered <- sf_summary %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))%>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

#basic filter 
sf_filtered_40 <- sf_summary %>%
  filter(trimmed.mean.quality > 40) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()



#filtering per length (L1)
sf_filtered_L1 <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & raw.length > 400) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()


#filtering per length (L2)
sf_filtered_L2 <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & raw.length > 400 & raw.length < 600) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()


#plot histogram of quality
p <- ggplot(data = sf_summary, aes(x = trim.start)) +
  geom_histogram(binwidth = 10)+
  scale_x_continuous(breaks = seq(0,600,50))+
  ylim(0,3000)

p1 <- ggplot(data = sf_summary, aes(x = trim.finish)) +
  geom_histogram(binwidth = 10) +
  scale_x_continuous(breaks = seq(0,600,50))+
  ylim(0,3000)

grid.arrange(p,p1)

#contiguos quality (trim_start & trim_finish)
sf_filtered_T <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & trim.start < 65 & trim.finish > 400) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

sf_summary %>% ggplot(aes(x = raw.mean.quality, y = raw.secondary.peaks)) +geom_point()
hist(sf_filtered_40$trim.finish)
#secondary peak counts

sf_filtered_S1 <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & trimmed.secondary.peaks < 50) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

sf_filtered_S2 <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & trimmed.secondary.peaks < 100) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()


#CDR3 tests

#here is the code used to create the listed S4 objects for secondary peaks detection
#do that for every sf_filtered (old filtered)
sf_filtered <- sf_40$summaries %>%
  filter(raw.length > 300 & raw.length < 600) %>%
  filter(trim.start < 65 & trim.finish > 150) %>%
  filter(trimmed.mean.quality > 30 & trimmed.secondary.peaks < 150) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*|\\_.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

pathnames_40 <- as.character(sf_filtered_40$file.path)
sangerseqlisted_40 <- sapply(pathnames, readsangerseq)
sp <- lapply(sangerseqlisted, secondary.peaks, ratio = 0.33)
df <- lapply(sp, function (x) x[["secondary.peaks"]])
df <- lapply(df, function(x) filter(x, position > 65 & position < 150))
df <- lapply(df, function(x) ifelse(nrow(x) > 0 , nrow(x), 0))
df <- tibble( sec.peak.CDR3 = as.numeric(as.character(df)), file.path = names(df))
df <- df %>% filter(sec.peak.CDR3 <= 5)

#CDR3 filter for summary
pathnames_summary <- as.character(sf_summary$file.path)
sangerseqlisted_summary <- sapply(pathnames_summary, readsangerseq)
sp_summary <- lapply(sangerseqlisted, secondary.peaks, ratio = 0.33)
df_summary<- lapply(sp, function (x) x[["secondary.peaks"]])
df_summary <- lapply(df, function(x) filter(x, position > 65 & position < 150))
df_summary <- lapply(df, function(x) ifelse(nrow(x) > 0 , nrow(x), 0))
df_summary <- tibble( sec.peak.CDR3 = as.numeric(as.character(df)), file.path = names(df))
df_summary <- df %>% filter(sec.peak.CDR3 <= 5)

#CDR3 filter for sf_filtered_40
pathnames_40 <- as.character(sf_filtered_40$file.path)
sangerseqlisted_40 <- sapply(pathnames_40, readsangerseq)
sp_40 <- lapply(sangerseqlisted, secondary.peaks, ratio = 0.33)
df_40 <- lapply(sp, function (x) x[["secondary.peaks"]])
df_40 <- lapply(df, function(x) filter(x, position > 65 & position < 150))
df_40 <- lapply(df, function(x) ifelse(nrow(x) > 0 , nrow(x), 0))
df_40 <- tibble( sec.peak.CDR3 = as.numeric(as.character(df)), file.path = names(df))
df_40 <- df %>% filter(sec.peak.CDR3 <= 5)

#plot histogram to check distribution of sec.peaks in CDR3
ggplot(df) +
  aes(x = sec.peak.CDR3) +
  geom_histogram(binwidth = 1) + 
  geom_text(stat = 'count', aes(label = stat(count), vjust = -0.2))

#creating final df for plottin
data.frame("initial" = 6282, 
  initial_no_repeat = 5664, 
  before = 3326, 
  before_65-150 = 2991, 
  filtered_40 = 4497,
  L1 = 4406,
  L2 = 4392,
  S1 = 4469,
  S2 = 4496,
  T = 3122)
col <- c("init", "init_NR", "P", "P_T", "filt_40", 
         "L1","L2","S1","S2","T", "C1", "C2")
val <- c(6282, 5664,3326,2991,4497,4406,4392,4469,4496,3122,3167,2057)

seq_n <- data.frame(col, val)
seq_n %>% arrange(desc(val)) %>% 
  mutate(col=factor(col, levels=col)) %>%
  ggplot(aes(x = col, y = val, label = val)) +
  geom_point() + geom_text_repel() +
  xlab("")+
  ylab("Number of filtered seqs")+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 60, vjust = .7, hjust = 1))
  




