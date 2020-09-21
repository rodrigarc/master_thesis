#program to create heatmap for different parameters using trim.cutoff as 20 
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