library(dplyr)
str(bm)

sars_like <- t(bm[c(11,12,13,19,21,22),])  

sars_like_0 <-  apply(row(sars_like)* !sars_like, 2, function(x) c(x[x!=0], x[x==0]))
sars_like_0_order <- sars_like_0[1:322,]
write.csv(sars_like_0_order, "sars_like_aminoacid_positions(0).csv")

hcov <- t(bm[c(14,15),]) 
hcov_1 <- as.data.frame(row(hcov) * hcov)

hcov_1_filter <- hcov_1 %>% filter(HCoV_OC43 == HCoV_HKU1 & HCoV_HKU1 > 0 & HCoV_OC43 > 0)
write.csv(hcov_1_filter, "hcov_aminoacid_positions(1).csv")

          