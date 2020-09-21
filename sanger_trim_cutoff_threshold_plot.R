#heatmap for trimmed_quality_score and trimmed_cutoff score
library(ggplot2)
library(hexbin)
library(cowplot)
library(dplyr)
library(RColorBrewer)

View(sf_20$summaries)
#data processing, combining dataframes etc.
sf20 <- sf_20$summaries %>% 
  mutate(sf = 20, 
         trimmed.integer = as.integer(trimmed.mean.quality)) %>%
  select(raw.mean.quality, trimmed.integer, sf, trimmed.mean.quality, trimmed.length)
sf20[is.na(sf20)] <- 0

sf30 <- sf_30$summaries %>% 
  mutate(sf = 30, 
         trimmed.integer = as.integer(trimmed.mean.quality)) %>%
  select(raw.mean.quality, trimmed.integer, sf, trimmed.mean.quality, trimmed.length)
sf30[is.na(sf30)] <- 0

sf40 <- sf_40$summaries %>% 
  mutate(sf = 40, 
         trimmed.integer = as.integer(trimmed.mean.quality)) %>%
  select(raw.mean.quality, trimmed.integer, sf, trimmed.mean.quality, trimmed.length)
sf40[is.na(sf40)] <- 0

sf50 <- sf_50$summaries %>% 
  mutate(sf = 50, 
         trimmed.integer = as.integer(trimmed.mean.quality)) %>%
  select(raw.mean.quality, trimmed.integer, sf, trimmed.mean.quality, trimmed.length)
sf50[is.na(sf50)] <- 0

combined <- rbind(sf20,sf30,sf40,sf50)

ggplot(combined, aes(x = sf, y = trimmed.integer, fill = raw.mean.quality))+
  geom_tile() + theme_cowplot()+
  xlab("Trim cutoff score")+
  ylab("Trimmed mean quality (integer)")

combined %>% group_by(sf) %>% 
  filter(trimmed.mean.quality > 50) %>% summarise(nrow = nrow,
                                                  filtering = 50,
                                                  raw_mean = mean(raw.mean.quality),
                                                  trimmed_mean = mean(trimmed.mean.quality),
                                                  mean_length = mean(trimmed_length),
                                                  sf)
# loop to count how many sequences were filtered using different trimmed.mean.quality
# trimmed.mean.quality from 1 to 60

datalist <- list()

qual_scores <- c(1:60)

for (i in qual_scores) {
  dat <- combined %>% 
    filter(trimmed.mean.quality > i) 
  dat$filtering <- i 
  dat <- dat %>% group_by(sf, filtering) %>%
    summarise(filtered_seq = length(filtering[filtering == i]),
              filtering = mean(filtering),
              raw_mean = mean(raw.mean.quality),
              trimmed_mean = mean(trimmed.mean.quality),
              mean_length = mean(trimmed.length),
              sf = mean(sf))
  datalist[[i]] <- dat # add it to your list
}

big_data <- do.call(rbind, datalist)

#trying graphs with filtered sequences
g1 <- ggplot(big_data, aes(x = sf, y = filtering, fill = raw_mean))+
  geom_tile() + theme_cowplot()+
  scale_fill_viridis_c()+
  xlab("Trim cutoff score")+
  ylab("Filtering score ( > x)")
g1
g2 <- ggplot(big_data, aes(x = sf, y = filtering, fill = filtered_seq))+
  geom_tile() + theme_cowplot()+
  scale_fill_viridis_c(option = "inferno")+
  xlab("Trim cutoff score")+
  ylab("Filtering score ( > x)")
g2

g3 <- ggplot(big_data, aes(x = sf, y = filtering, fill = as.integer(trimmed_mean)))+
  geom_tile() + theme_cowplot()+
  xlab("Trim cutoff score")+
  ylab("Filtering score ( > x)")+
  labs(fill = "Trimmed_mean")+
  scale_fill_viridis_c(option = "magma")+
  ylim(0,60)

g4 <- ggplot(big_data, aes(x = sf, y = filtering, fill = mean_length))+
  geom_tile() + theme_cowplot()+
  xlab("Trim cutoff score")+
  ylab("Filtering score ( > x)")+
  labs(fill = "mean length")+
  scale_fill_viridis_c(option = "magma")+
  ylim(0,60)
g4
gridExtra::grid.arrange(g1,g2, g3, g4)

ggplot(big_data, aes(x = sf, y = trimmed_mean, fill = raw_mean)) +
  geom_hex(bins = 50) +
  scale_fill_continuous(type = "viridis") +
  theme_cowplot()

ggplot(big_data, aes(x = sf, y = trimmed_mean, fill = raw.mean)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white")+
  theme_cowplot()

ggplot(big_data, aes(x = sf, y = trimmed_mean, fill = raw_mean)) + geom_point()

