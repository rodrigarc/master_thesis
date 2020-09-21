library(dplyr)
library(data.table)
library(sangerseqR)
library(stringr)
library(Biostrings)

getwd()
setwd("~/Box Sync")




#load sf_trim_0.01 (10-20 trim)

sf_20 <- sf_trim_0.01$summaries %>% 
  mutate(sf = 20, 
         trimmed.integer = as.integer(trimmed.mean.quality)) 
sf_20[is.na(sf_20)] <- 0

#naming
sf_trim_0.01$summaries <- sf_trim_0.01$summaries %>%
  mutate(plate = paste(substring(folder.name,5,6)))%>%
  mutate(well = gsub("\\-.*|\\_.*", "", file.name)) %>%
  mutate(ID = substring(folder.name, 1, 3)) %>%
  mutate(chain = "HC")%>% 
  mutate(file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_")))
sf_trim_0.01$summaries[is.na(sf_trim_0.01$summaries)] <- 0

#no filter
sf_full_dedup <- sf_trim_0.01$summaries %>%
  group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

#old filter
sf_filtered_old_dedup <- sf_trim_0.01$summaries %>%
    filter(raw.length > 300 & raw.length < 600) %>%
    filter(trim.start < 65 & trim.finish > 150) %>%
    filter(trimmed.mean.quality > 30 & trimmed.secondary.peaks < 150) %>%
    mutate(file_ID = sub("_R","",paste(ID,plate,well,chain,sep = "_"))) %>%
    group_by(file_ID) %>%
    filter(raw.mean.quality == max(raw.mean.quality)) %>%
  ungroup()

#testing filters
sf_filtered_4f <- sf_trim_0.01$summaries %>%
  filter(raw.length >= 400 &
           trim.start < 33 &
           trim.finish > 409 &
           trimmed.mean.quality >= 30 )

sf_filtered_4f_peak <- sf_trim_0.01$summaries %>%
  filter(raw.length >= 400 &
           trim.start < 33 &
           trim.finish > 409 &
           trimmed.mean.quality >= 30 &
           trimmed.secondary.peaks <= 50 )


#adding cdr3 filter

sf_filtered_4f_peak_cdr3 <- sf_trim_0.01$summaries %>%
  filter(raw.length >= 400 &
           trim.start < 33 &
           trim.finish > 409 &
           trimmed.mean.quality >= 30 &
           trimmed.secondary.peaks <= 50 )

#pathnames <- as.character(sf_full_dedup$file.path)
#sangerseqlisted <- sapply(pathnames, readsangerseq)
#sp <- lapply(sangerseqlisted, secondary.peaks, ratio = 0.33)
#df <- lapply(sp, function (x) x[["secondary.peaks"]])
#df <- lapply(df, function(x) filter(x, position > 100 & position < 150))
#df <- lapply(df, function(x) ifelse(nrow(x) > 0 , nrow(x), 0))
#df <- tibble( sec.peak.CDR3 = as.numeric(as.character(df)), file.path = names(df))
#df <- df %>% filter(sec.peak.CDR3 <= 5)

#read df with cdr3
df <- read.csv("~/Box Sync/codes/R filtering improvement/CDR3_sec_pek_df.csv")
sf_dedup_cdr3 <- merge(sf_dedup, df)

#sf filtered and deduplicated
3882-sf_filtered_4f_peak_cdr3 %>% group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>% nrow()

3882- sf_filtered_4f %>% group_by(file_ID) %>%
  filter(raw.mean.quality == max(raw.mean.quality)) %>% nrow()

sf_filtered_4f_dedup <- sf_filtered_4f %>% group_by(file_ID) %>%
  filter(trimmed.mean.quality == max(trimmed.mean.quality)) %>%
  ungroup()
sf_filtered_4f_peak_dedup <- sf_filtered_4f_peak %>% group_by(file_ID) %>%
  filter(trimmed.mean.quality == max(trimmed.mean.quality)) %>%
  ungroup()
sf_filtered_4f_peak_cdr3_dedup <- sf_filtered_4f_peak_cdr3 %>% group_by(file_ID) %>%
  filter(trimmed.mean.quality == max(trimmed.mean.quality)) %>%
  ungroup()
sf_filtered_old_dedup <- sf_filtered_old %>% group_by(file_ID) %>%
  filter(trimmed.mean.quality == max(trimmed.mean.quality)) %>%
  ungroup()


#create fasta files
filtered.filepath <- as.character(sf_full_dedup$file.path)
sangerseqfiltered <- sapply(filtered.filepath, readsangerseq)
sangerbasecallfasta <- sapply(sangerseqfiltered, primarySeq)
#name_fasta <- sf_filtered_4f %>% mutate(name = file_ID) 
name_fasta <- sf_full_dedup$file_ID
names(sangerbasecallfasta) <- name_fasta

ape::write.dna(sangerbasecallfasta, 
               file = "~/Box Sync/RSV NGS/sanger_QC/fasta-filters/sanger_QC-full_dedup.fasta", 
               format = 'fasta', 
               nbcol = -1, 
               colsep = "",
               colw = 10000000,
               append = F)


#not using phread score anymore
cdr3 <- lapply(sf_trim_0.01$quality_scores, function(x)
  filter(x,x[["position"]] >= 70 & x[["position"]] <= 160 & x[["score"]] <= 5))

delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
  x.list[unlist(lapply(x.list, length) != 0)]
}
cdr3_df <- rbindlist(cdr3, idcol = T)
cdr3_df <- cdr3_df %>% rename(file.path = .id)
unique(cdr3_df$file.path)
str(sf_filtered_4f_peak[!sf_filtered_4f_peak$file.path %in% cdr3_df$file.path,])

sf_filtered_4f_peak_cdr3 <- sf_filtered_4f_peak[!(sf_filtered_4f_peak$file.path %in% cdr3_df$file.path),]

delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
  x.list[unlist(lapply(x.list, length) != 0)]
}



#plotting cdr3
library(ggplot2)
cdr3 <-  sf_trim_0.01$quality_scores
cdr3_b <- rbindlist(cdr3)

g <- cdr3_b %>% ggplot(aes(x = position, y = score)) 
g1 <- g + geom_smooth(method = "gam")
g1
g2 <- g1 + geom_point()
g2

cdr3_s <- cdr3_b %>% group_by(position) %>% summarise(scores = mean(score),
                                                      stand = sd(score))
cdr3_s %>% filter(position >= 50 & position <=150) %>%
  ggplot(aes(position, scores)) + geom_line() + 
  geom_ribbon(aes(ymin = scores+stand, 
                  ymax = scores-stand, 
                  fill = "red"), 
              alpha = .15) + theme_bw()
cdr3_s %>% filter(position >= 50 & position <=150) %>%
  ggplot(aes(position, scores)) + geom_line()


#reading filtered tabs and checking the effect of each filter
ls <- list.files("~/Box Sync/RSV NGS/sanger_QC/igdiscover/filtered-tabs", 
                 full.names = T)
ls <- grep("dedup", ls, value = T)
tabs <- lapply(ls, fread)
tabs <- lapply(tabs, as.data.frame)
names(tabs) <- ls

#4f - 4f_peak
a <- dplyr::setdiff(tabs[[1]][["name"]],
                    tabs[[3]][["name"]])

#4f - 4f_peak_cdr3
b <- dplyr::setdiff(tabs[[1]][["name"]],
                    tabs[[2]][["name"]])

#4_peak - 4f_peak_cdr3
c <- dplyr::setdiff(tabs[[3]][["name"]],
                    tabs[[2]][["name"]])
#merging
a_score <- subset(sf_full_dedup, file_ID %in% a)
b_score <- subset(sf_full_dedup, file_ID %in% b)
c_score <- subset(sf_full_dedup, file_ID %in% c)

merged <- list(a_score,b_score,c_score)
names(merged) <- c("a","b","c")
merged <- rbindlist(merged, idcol = T)

pathnames <- as.character(merged$file.path)
sangerseqlisted <- sapply(pathnames, readsangerseq)
sp <- lapply(sangerseqlisted, secondary.peaks, ratio = 0.33)
df <- lapply(sp, function (x) x[["secondary.peaks"]])
df <- lapply(df, function(x) filter(x, position > 100 & position < 150))
df <- lapply(df, function(x) ifelse(nrow(x) > 0 , nrow(x), 0))
df <- tibble( sec.peak.CDR3 = as.numeric(as.character(df)), file.path = names(df))
df <- df %>% filter(sec.peak.CDR3 > 5)

merged_t <- merge(merged, df)
merged_t <- merged_t %>% group_by(file_ID)

merged_summary <- merged_t %>% group_by(.id) %>% summarise(mean_raw_length = mean(raw.length), 
                                       sd_raw_length = sd(raw.length),
                                       mean_raw_quality = mean(raw.mean.quality), 
                                       sd_raw_quality = sd(raw.mean.quality),
                                       mean_trimmed_quality = mean(trimmed.mean.quality), 
                                       sd_trimmed_quality = sd(trimmed.mean.quality),
                                       mean_secondary_peaks = mean(trimmed.secondary.peaks), 
                                       sd_secondary_peaks = sd(trimmed.secondary.peaks),
                                       mean_cdr3_peaks = mean(sec.peak.CDR3 ),
                                       sd_cdr3_peaks = sd(sec.peak.CDR3),
                                       mean_trim_start = mean(trim.start), 
                                       sd_trim_start = sd(trim.start),
                                       mean_trim_finish = mean(trim.finish), 
                                       sd_trim_finish = sd(trim.finish))
write.csv(merged_summary, 
          "~/Box Sync/RSV NGS/sanger_QC/data/2020-09-08/summary_filtering_scores.csv")
