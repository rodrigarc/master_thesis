library(stringr)
library(dplyr)

tb <- Diversity_all_table %>% 
  mutate(animal = gsub("\\_.*", "", sample_name),
         spec = gsub("*._", "", sample_name))


tb <- Diversity_all_table %>% 
  mutate(animal = gsub("\\_.*", "", sample_name),
         spec = ifelse(grepl("total", sample_name), "Total", 
                       ifelse(grepl("pref",sample_name), "PreF", 
                              ifelse(grepl("dp", sample_name), "DP", ""))))
