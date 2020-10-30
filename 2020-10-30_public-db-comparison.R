library(tidyverse)

test <- openxlsx::read.xlsx("~/Box Sync/RSV NGS/201030 Database comparison.xlsx")

test_omit <- na.omit(test[2:11])

test_reshaped <- test_omit %>% gather() 

ggplot(gather(test_reshaped), aes(value)) + 
  geom_histogram(bins = 10) + 
  facet_wrap(~key, scales = 'free_x')
