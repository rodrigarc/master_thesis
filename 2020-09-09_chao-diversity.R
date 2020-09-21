library(iNEXT)
library(vegan)
library(dplyr)
library(reprex)
library(ggplot2)
library(data.table)
library(reshape2)
library(cowplot)
require(openxlsx)

getwd()
setwd("~/Desktop/sebastian/replicates-NGS/")

#read the files 
#####check if it is read.csv or read.csv2, it was confusing before
pref <- as.matrix(read.csv2("~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clones 4 Vegan+Recon/Clones_PreF.csv", header = F))
dp <- as.matrix(read.csv2("~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clones 4 Vegan+Recon/Clones_DP.csv",  header = F))
total <- as.matrix(read.csv("~/Desktop/rodrigo/vegan/V4_NGS_TR1_Clones 4 Vegan/Clones_total.csv",  header = F))

#change NA to 0
pref[is.na(pref)] <- 0
dp[is.na(dp)] <- 0
total[is.na(total)] <- 0


#create function to subsample based on y value, replicate 100 times, 
#calculate the estimators, then the mean of the replicates
chaox100 <- function(x){
  replicate(100, {
  subsample <-  rrarefy(x,min(rowSums(x)))
  chao <- estimateR(subsample)
})}

#using it, y it is the minimum number of clone, I just check on excel the animal
#with the lowest number of columns and checked how many columns it had
chao_pref <- chaox100(pref)
chao_dp <- chaox100(dp)
chao_total <- chaox100(total)

chao_pref_list <- lapply(asplit(chao_pref,1),t)
chao_dp_list <- lapply(asplit(chao_dp,1),t)
chao_total_list <- lapply(asplit(chao_total,1),t)


chao_pref_df <- bind_rows(chao_pref_list)
chao_dp_df <- do.call(rbind, chao_dp_list)
chao_total_df <- data.frame(do.call(rbind, chao_total_list))

s.obs_dp  <- as.data.frame((chao_dp_list$S.obs))
s.obs_dp$id<- rownames(s.obs_dp)
g3 <- ggplot(data = melt(s.obs_dp), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. observed")+
  xlab("ID") +
  theme(legend.position = "none") +
  ylim(0,200)+
  ggtitle("DP")
g3
s.obs_total  <- as.data.frame((chao_total_list$S.obs))
s.obs_total$id<- rownames(s.obs_total)
g2 <- ggplot(data = melt(s.obs_total), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. observed")+
  xlab("ID") +
  theme(legend.position = "none") +
  ylim(0,200)+
  ggtitle("Total")
g2
s.obs_pref  <- as.data.frame((chao_pref_list$S.obs))
s.obs_pref$id<- rownames(s.obs_pref)
g1 <- ggplot(data = melt(s.obs_pref), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. observed")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("PreF")
g1 

s.chao_pref  <- as.data.frame((chao_pref_list$S.chao1))
s.chao_pref$id<- rownames(s.chao_pref)
g4 <- ggplot(data = melt(s.chao_pref), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. chao1")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("PreF")
g4

s.chao_total  <- as.data.frame((chao_total_list$S.chao1))
s.chao_total$id<- rownames(s.chao_total)
g5 <- ggplot(data = melt(s.chao_total), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. chao1")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("Total")

s.chao_dp  <- as.data.frame((chao_dp_list$S.chao1))
s.chao_dp$id<- rownames(s.chao_dp)
g6 <- ggplot(data = melt(s.chao_dp), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. chao1")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("DP")

s.ace_dp  <- as.data.frame((chao_dp_list$S.ACE))
s.ace_dp$id<- rownames(s.ace_dp)
g7 <- ggplot(data = melt(s.ace_dp), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. ACE")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("DP")

s.ace_pref  <- as.data.frame((chao_pref_list$S.ACE))
s.ace_pref$id<- rownames(s.ace_pref)
g8 <- ggplot(data = melt(s.ace_pref), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. ACE")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("PreF")

s.ace_total  <- as.data.frame((chao_total_list$S.ACE))
s.ace_total$id<- rownames(s.ace_total)
g9 <- ggplot(data = melt(s.ace_total), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) +
  scale_fill_viridis_d()+
  theme_cowplot()+
  ylab("S. ACE")+
  xlab("ID") +
  theme(legend.position = "none")+
  ylim(0,200)+
  ggtitle("Total")


gridExtra::grid.arrange(g2,g1,g3, g5, g4, g6, g9,g8,g7, ncol = 3, nrow = 3)


#calculate means
list_names <- c("s.ace_dp", 
                "s.ace_pref", 
                "s.ace_total", 
                "s.chao_dp",
                "s.ace_pref",
                "s.chao_total",
                "s.obs_dp", 
                "s.obs_pref", 
                "s.obs_total")
list_estimators <- list(s.ace_dp, 
                        s.ace_pref, 
                        s.ace_total, 
                        s.chao_dp,
                        s.chao_pref,
                        s.chao_total,
                        s.obs_dp, 
                        s.obs_pref, 
                        s.obs_total)
names(list_estimators) <- list_names

#save in excel worksheets
list_of_datasets <- list("s.ace_dp" = s.ace_dp, 
                         "s.ace_pref" = s.ace_pref, 
                         "s.ace_total" = s.ace_total, 
                         "s.chao_dp" = s.chao_dp,
                         "s.chao_pref" = s.ace_pref,
                         "s.chao_total" = s.chao_total,
                         "s.obs_dp" = s.obs_dp, 
                         "s.obs_pref" = s.obs_pref, 
                         "s.obs_total" = s.obs_total)
#calculate mean in excel to each estimator per animal
write.xlsx(list_of_datasets, 
           file = "~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clones 4 Vegan+Recon/Vegan output/chao_datasets.xlsx")
