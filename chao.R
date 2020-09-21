library(vegan)
library(SpadeR)
library(dplyr)

#read the files
pref <- as.matrix(read.csv("~/Downloads/scClonotypes4vegan-PreF.csv", header = F))
dp <- as.matrix(read.csv("~/Downloads/scClonotypes4vegan-DP.csv",  header = F))

#change NA to 0
pref[is.na(pref)] <- 0
dp[is.na(dp)] <- 0

#create function to subsample based on y value, replicate 1000 times, 
#calculate the estimators, then the mean of the replicates
chaox100 <- function(x,y){
  replicate(1000, {
  subsample <-  rrarefy(x,y)
  chao <- ChaoSpecies(subsample)
  df <- chao$Species_table
}) %>% apply(MARGIN = 1:2, mean)}

#using it, y it is the minimum number of clone, I just check on excel the animal
#with the lowest number of columns and checked how many columns it had
chao_pref <- chaox100(pref,70)
chao_dp <- chaox100(dp, 64)

write.csv(chao_pref,"~/Desktop/rodrigo/chao/chao_pref.csv")
write.csv(chao_dp,"~/Desktop/rodrigo/chao/chao_dp.csv")
