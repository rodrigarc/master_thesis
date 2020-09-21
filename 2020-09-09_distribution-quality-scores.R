library(dplyr)

#naming
sf_trim_0.01$summaries <- sf_trim_0.01$summaries %>%
  mutate(plate = paste(substring(folder.name,5,6)))%>%
  mutate(well = gsub("\\-.*|\\_.*", "", file.name)) %>%
  mutate(ID = substring(folder.name, 1, 3)) %>%
  mutate(chain = "HC")%>% 
  mutate(file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))
sf_trim_0.01$summaries[is.na(sf_trim_0.01$summaries)] <- 0

#dedup
sf_dedup <- sf_trim_0.01$summaries %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()
#read df with cdr3
df <- read.csv("~/Box Sync/codes/R filtering improvement/CDR3_sec_pek_df.csv")
sf_dedup_cdr3 <- merge(sf_dedup,df)

length <- as.data.frame(sf_dedup_cdr3$raw.length)
sf_dedup_l <- sf_dedup_cdr3 %>% filter(length >= 400) %>%
  select(sec.peak.CDR3,trim.start, trim.finish, trimmed.secondary.peaks, trimmed.mean.quality)

write.csv(sf_dedup_l, "~/Box Sync/RSV NGS/sanger_QC/data/2020-09-09/distribution_quality_scores.csv",
          row.names = F)
write.csv(sf_dedup_l, "~/Box Sync/Gunilla presentation/distribution_quality_scores.csv",
          row.names = F)
write.csv(length, "~/Box Sync/RSV NGS/sanger_QC/data/2020-09-09/distribution_length.csv",
          row.names = F)
write.csv(length, "~/Box Sync/Gunilla presentation/distribution_length.csv",
          row.names = F)



