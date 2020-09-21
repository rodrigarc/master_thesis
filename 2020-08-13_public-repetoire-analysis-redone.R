#rename single-cell clontypes by number
#query using clonotypes file (which is filtered tab)



library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#after rerunning clonoquery with a unique name to each clone
ls <- list.files("~/Desktop/rodrigo/public-clonotype/Public_clonoquery_3/IgM/", 
                full.names = T)
ls <- grep("summary", ls, value = T)
ls
x <- lapply(ls, fread, fill = T)
names(x) <- c("E11", "E12", "E14","E16","E17", "E18", "E21", "E23","E24")

test_3 <- rbindlist(x, idcol = TRUE)
test_3_s <- test_3 %>% select(-c(18,19))
test_3_f <- test_3_s %>% filter(size > 0 ) %>%
  group_by(name) %>% 
  summarise(size_rep_seq = sum(size), animals_rep_seq = paste0(unique(.id), collapse = ", ")) 
test_n <- test_3_f %>%  
  mutate(number_animals_rep_seq = str_count(animals_rep_seq, "E")) %>%
  filter(number_animals_rep_seq > 1) %>%
  arrange(desc(number_animals_rep_seq))


#using clonoquery full to get all the sequences names

y <- fread("~/Desktop/rodrigo/public-clonotype/Public_clonotypes/2020-08-10_public_clonotypes_full.txt", fill = T)
df <- y %>%
  mutate(clone = paste0(cumsum(name == "") + 1),
         clone = if_else(name == "", "", clone)) %>%
  mutate(ID = sub("\\_.*","", name)) 
df <- df[!apply(is.na(df) | df == "", 1, all),]

df <- df %>% rename(sequence = name, 
                    name = clone,
                    size = count)
df$name <- paste0("P",sprintf("%04s", df$name))

df_g <- df %>% group_by(name) %>% 
  summarise(size_sc = sum(size),
            animals_sc = paste0(unique(ID), collapse = ", "),
            sequence = paste0(sequence, collapse = ", "))
df_sc <- df_g %>%  
  mutate(number_animals_sc = str_count(df_g$animals_sc, "E")) %>%
  arrange(desc(number_animals_sc))
merged <- merge(df_sc,test_n, by = "name") 
merged_a <- merged %>% arrange(desc(number_animals_rep_seq)) %>% 
  rename(clone = name) %>%
  select(clone, animals_rep_seq, animals_sc, 
         number_animals_rep_seq, number_animals_sc,
         size_rep_seq, size_sc, sequence)
  


#csv files from dfs
write.csv(merged_a, 
          "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v3/2020-08-14_public_table_IgM_summary.csv",
          row.names = F)


library(kableExtra)
library(knitr)
knitr::kable(merged_a) %>%
  kable_styling(bootstrap_options = "striped", full_width = F)




