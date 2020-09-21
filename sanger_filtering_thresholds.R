library(dplyr)
library(gridExtra)

#sanger_filtering_thresholds_tests

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
         well = gsub("\\-.*", "", file.name),
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
         well = gsub("\\-.*", "", file.name),
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
         well = gsub("\\-.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))%>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

#plot histogram of quality
p <- ggplot(data = sf_summary, aes(x = trim.start)) +
  geom_histogram(binwidth = 10)+
  scale_x_continuous(breaks = seq(0,600,50))+
  ylim(0,3000)+
  t

p1 <- ggplot(data = sf_summary, aes(x = trim.finish)) +
  geom_histogram(binwidth = 10) +
  scale_x_continuous(breaks = seq(0,600,50))+
  ylim(0,3000)

grid.arrange(p,p1)

#contiguos quality (trim_start & trim_finish)
sf_filtered_T <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & trim.start < 65 & trim.finish > 400) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))%>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

#secondary peak counts

sf_filtered_S1 <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & trimmed.secondary.peaks < 50) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))%>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

sf_filtered_S2 <- sf_summary %>%
  filter(trimmed.mean.quality > 40 & trimmed.secondary.peaks < 100) %>%
  mutate(plate = paste(substring(folder.name,5,6)),
         well = gsub("\\-.*", "", file.name),
         ID = substring(folder.name, 1, 3),
         chain = "HC",
         file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))%>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()
