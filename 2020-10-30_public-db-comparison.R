library(tidyverse)
library(data.table)
library(rstatix)
library(ggpubr)


ls <- list.files("~/Box Sync/RSV NGS/Database comparison/filtered_tabs/", full.names = T)
ls_names <- list.files("~/Box Sync/RSV NGS/Database comparison/filtered_tabs/")

tabs <- lapply(ls, fread)
tabs_v <- lapply(tabs, select, V_SHM, name)
names(tabs_v) <- ls_names
tabs_v_s <- rbindlist(tabs_v,idcol = T)


####### checking histograms, qq plots and test for normality using shapiro_test ######
tabs_v_s %>% ggplot(aes(V_SHM)) +
  geom_histogram(bins = 10) + 
  facet_wrap(~.id)
ggplot(tabs_v_s, aes(sample=V_SHM)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap( ~ .id)
tabs_v_s %>% group_by(.id) %>% shapiro_test(V_SHM)


tabs_v_s <- tabs_v_s %>% mutate(pdb = ifelse(grepl("individualized", .id), "individualized", "public"))

stat.test <- tabs_v_s %>% t_test(V_SHM ~ .id, p.adjust.method = "holm") 

tabs_v_s$.id <- sub("-filtered.tab", "",x = tabs_v_s$.id)

tabs_v_s$.id <- factor(tabs_v_s$.id, levels = c( "IMGThm","IMGTrm", 
                                                 "ramesh", "Cirelli", 
                                                 "sundling", "Corcoran", 
                                                 "Francica", "publicdb", "mmul10",  "zhang",
                                                "individualized.tab"))
tabs_v_s %>% group_by(.id) %>% summarise(avg = mean(V_SHM)) %>% arrange(avg)



# Create the plot

myplot <- ggviolin(
  tabs_v_s, x = ".id", y = "V_SHM",
  fill = ".id", palette = "npg", legend = "none",
  ggtheme = theme_pubr(border = TRUE)
) + labs(y = "SHM %", x = "")  + rotate_x_text(angle = 65)
myplot

# Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = ".id") %>% filter(group1 == "individualized.tab"|group2 == "individualized.tab") %>% 
  as.data.frame()
class(stat.test)
#write.csv(stat.test[-c(11:14)], "~/Box Sync/RSV NGS/Database comparison/plot/2020-02-20_t-test_p-adjusted.csv",
          row.names = F)


myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif", bracket.nudge.y = -15) +
  rotate_x_text(angle = 65)



##create a new public_db

