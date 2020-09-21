library(stringr)
library(dplyr)
library(openxlsx)

#add this script in terminal to create diversity table
$ python2 recon_v2.5.py -D -b error_bar_parameters.txt -o Diversity_total_table.txt 1_total_fitfile.txt 2_total_fitfile.txt 3_total_fitfile.txt 4_total_fitfile.txt 5_total_fitfile.txt 6_total_fitfile.txt 7_total_fitfile.txt 8_total_fitfile.txt 9_total_fitfile.txt 1_dp_fitfile.txt 2_dp_fitfile.txt 3_dp_fitfile.txt 4_dp_fitfile.txt 5_dp_fitfile.txt 6_dp_fitfile.txt 7_dp_fitfile.txt 8_dp_fitfile.txt 9_dp_fitfile.txt 1_pref_fitfile.txt 2_pref_fitfile.txt 3_pref_fitfile.txt 4_pref_fitfile.txt 5_pref_fitfile.txt 6_pref_fitfile.txt 7_pref_fitfile.txt 8_pref_fitfile.txt 9_pref_fitfile.txt


Diversity_all_table <- read.table("~/Desktop/Recon/2020-09-09/Diversity_total_table.txt", header = T)

tb <- Diversity_all_table %>% 
  mutate(animal = gsub("\\_.*", "", sample_name),
         spec = gsub("*._", "", sample_name))


tb <- Diversity_all_table %>% 
  mutate(animal = gsub("\\_.*", "", sample_name),
         spec = ifelse(grepl("total", sample_name), "Total", 
                       ifelse(grepl("pref",sample_name), "PreF", 
                              ifelse(grepl("dp", sample_name), "DP", "")))) %>%
  relocate(animal, spec)
write.xlsx(tb, 
           "~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clones 4 Vegan+Recon/Recon output/200909_diversity_estimation_recon.xlsx")










