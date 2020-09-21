#iNEXT

library(iNEXT)
library(dplyr)
library(reprex)
library(ggplot2)
library(data.table)
library(reshape2)
library(cowplot)
library(openxlsx)
library(gridExtra)

getwd()
setwd("~/Desktop/sebastian/replicates-NGS/")

#read the files
pref <- as.matrix(read.csv("~/Desktop/rodrigo/vegan/V4_NGS_TR1_Clones 4 Vegan/Clones_PreF.csv", header = F))
dp <- as.matrix(read.csv("~/Desktop/rodrigo/vegan/V4_NGS_TR1_Clones 4 Vegan/Clones_DP.csv",  header = F))
total <- as.matrix(read.csv("~/Desktop/rodrigo/vegan/V4_NGS_TR1_Clones 4 Vegan/Clones_total.csv",  header = F))

#change NA to 0
pref[is.na(pref)] <- 0
dp[is.na(dp)] <- 0
total[is.na(total)] <- 0

#testing iNEXT
data(spider)
str(spider)

#PreF plot
list_pref <- split(pref, row(pref))
i_pref <- iNEXT(list_pref, nboot = 100)
g1 <- ggiNEXT(i_pref) + theme_cowplot()+ 
  scale_shape_manual(values=c(19,20,15,16,8,17,18,9,10))+
  ggtitle("PreF")
g1
#DP plot
list_dp <- split(dp, row(dp))
i_dp <- iNEXT(list_dp, nboot = 100)
g2 <- ggiNEXT(i_dp) + theme_cowplot()+ 
  scale_shape_manual(values=c(19,20,15,16,8,17,18,9,10))+
  ggtitle("DP")

#Total plot
list_total <- split(total, row(total))
i_total <- iNEXT(list_total, nboot = 100)
g3 <- ggiNEXT(i_dp) + theme_cowplot()+ 
  scale_shape_manual(values=c(19,20,15,16,8,17,18,9,10))+
  ggtitle("Total")
g3+legend

g1
g2
g3

gridExtra::grid.arrange(g3,g1,g2, nrow = 1)
grid.arrange(g3,g1,g2, nrow=1) + opts(legend.position="bottom")

#extracting values
write.xlsx(i_total$AsyEst, "iNext_total.xlsx")
write.xlsx(i_dp$AsyEst, "iNext_DP.xlsx")
write.xlsx(i_pref$AsyEst, "iNext_PreF.xlsx")


