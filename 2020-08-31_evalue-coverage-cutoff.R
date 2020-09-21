library(dplyr)
library(openxlsx)
library(data.table)
library(ggplot2)

#reading all files containing evalue and coverage to set the threshold to flag
#to avoid having to read it again, Rdata was saved in folder NGS RSV > mAb_QC
sc_hc <- read.xlsx(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_Single cell/sc_v4 filtered files/sc_summary_filtered.xlsx")
sc_lc <- read.xlsx("~/Box Sync/RSV NGS/v4_Analysis/LCs/200806_LC_all.xlsx")

ls_B1 <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/rep-seq_filtered-files/B1/", full.names = T)
rs_B1 <- lapply(ls_B1, fread, select = c("V_covered", "J_covered", "V_evalue", "J_evalue"))
rs_B1 <- rbindlist(rs_B1)

ls_B2 <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/rep-seq_filtered-files/B2/", full.names = T)
rs_B2 <- lapply(ls_B2, fread, select = c("V_covered", "J_covered", "V_evalue", "J_evalue"))
rs_B2 <- rbindlist(rs_B2)

ls_igm <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/rep-seq_filtered-files/IgM", full.names = T)
rs_igm <- lapply(ls_igm, fread, select = c("V_covered", "J_covered", "V_evalue", "J_evalue"))
rs_igm <- rbindlist(rs_igm)

combined <- rbindlist(list(B1 = rs_B1, 
                           B2 = rs_B2, 
                           IgM = rs_igm, 
                           sc_HC = sc_hc[c("V_covered", "J_covered", "V_evalue", "J_evalue")],
                           sc_LC = sc_lc[c("V_covered", "J_covered", "V_evalue", "J_evalue")]), 
                      idcol = T)
str(combined)

#plot boxplots from all the datasets
g <- combined %>% ggplot(aes(x = .id, y = log10(V_evalue))) + geom_boxplot() + xlab("")
g1 <- combined %>% ggplot(aes(x = .id, y = V_covered)) + geom_boxplot() + xlab("")
g2 <- combined %>% ggplot(aes(x = .id, y = log10(J_evalue))) + geom_boxplot() + xlab("")
g3 <- combined %>% ggplot(aes(x = .id, y = J_covered)) + geom_boxplot() + xlab("")

gridExtra::grid.arrange(g,g1,g2,g3, nrow = 2)

#calculate threshold to flag
str(combined)
statistic <- combined 
#testing the parameters for outlier
stats_threshold_outlier <- statistic %>% group_by(.id) %>% 
  summarise(V_covered_lower = -(1.5*IQR(V_covered)) + quantile(V_covered)[2],
            V_covered_upper = (1.5*IQR(V_covered)) + quantile(V_covered)[4],
            V_evalue_lower = -(1.5*IQR(V_evalue)) + quantile(V_evalue)[2],
            V_evalue_upper = (1.5*IQR(V_evalue)) + quantile(V_evalue)[4],
            J_covered_lower = -(1.5*IQR(J_covered)) + quantile(J_covered)[2],
            J_covered_upper = (1.5*IQR(J_covered)) + quantile(J_covered)[4],
            J_evalue_lower = -(1.5*IQR(J_evalue)) + quantile(J_evalue)[2],
            J_evalue_upper = (1.5*IQR(J_evalue)) + quantile(J_evalue)[4])

#testing the parameters for 75% quantile
stats_threshold_75th <- statistic %>% group_by(.id) %>% 
  summarise(V_covered_lower = quantile(V_covered)[2],
            V_evalue_upper = quantile(V_evalue)[4],
            J_covered_lower = quantile(J_covered)[2],
            J_evalue_upper = quantile(J_evalue)[4])

stats_threshold_75th
