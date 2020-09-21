#rename single-cell clontypes by number
#query using clonotypes file (which is filtered tab)



#rethink how to calculate the distance between each clonotype and 
#download clonoquery TR1_v4 for every animal and then calculate
#between every animal
#plot using igraph

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#after rerunning clonoquery with a unique name to each clone
ls <- list.files("~/Desktop/rodrigo/public-clonotype/Public_clonoquery_3/IgM", 
                full.names = T)
ls <- grep("summary", ls, value = T)
ls
x <- lapply(ls, fread, fill = T)
names(x) <- c("E11", "E12", "E14","E16","E17", "E18", "E21", "E23","E24")

test_3 <- rbindlist(x, idcol = TRUE)
test_3_s <- test_3 %>% select(-c(18,19))
test_3_f <- test_3_s %>% filter(size > 0 ) %>%
  group_by(name) %>% 
  summarise(size = sum(size), animals = paste0(unique(.id), collapse = ", ")) %>%
  filter(nchar(animals) > 4)
test_n <- test_3_f %>%  
  mutate(number_animals_rep_seq = str_count(test_3_f$animals, "E")) %>%
  arrange(desc(number_animals_rep_seq))


#combine clonoquery in the rep-seq with the sc data
test_3_f <- test_3_f %>% rename(size_rep_seq = size, 
                                animals_rep_seq = animals)
sc <- read.csv("~/Desktop/rodrigo/public-clonotype/Public_clonotypes/2020-08-10_public_clonotypes_query_ref-table.csv")
sc <- sc %>% rename(sequence = name, 
                    name = Public.clonotype,
                    size = count)
sc <- sc %>% mutate(ID = gsub("(.+?)(\\_.*)", "\\1", sequence))

sc <- sc %>%
  group_by(name) %>% 
  summarise(size_sc = sum(size),
            animals_sc = paste0(unique(ID), collapse = ", "),
            sequence)
merged <- merge(test_3_f, sc, by = "name")

#not needed to delete empty rows because here we are using summary
#x <- lapply(x, function (x) x[!apply(is.na(x) | x == "", 1, all),])       
x_all <- do.call(rbind, x)

#getting the clonoquery named properly, all clones queries, aka full table
x_m_full <- x_all %>% 
  mutate(grp = cumsum(grepl("Query",name))) %>%
  group_by(grp) %>% 
  mutate(clonoquery = ifelse(row_number() > 1, first(name), name)) %>%
  ungroup()

#filtering only to the clones that were found, aka a summary table
x_m_summary <- x_all %>% 
  mutate(grp = cumsum(grepl("Query",name))) %>%
  group_by(grp) %>% 
  mutate(clonoquery = ifelse(row_number() > 1, first(name), "")) %>%
  ungroup() %>%
  filter(clonoquery != "") %>%
  select(-c(grp)) %>%
  relocate(clonoquery)

#writing table of public clonotypes
x_m_basic_table <- x_m_summary %>%
  mutate(ID = name) %>%
  group_by(clonoquery) %>%
  summarise(animals = paste0(unique(ID), collapse = ", "))

library(kableExtra)
library(knitr)
knitr::kable(test_n) %>%
kable_styling(bootstrap_options = "striped", full_width = F)
  


#csv files from dfs
write.csv(x_all, 
          "~/Desktop/rodrigo/public-clonotype/Analysis/Public_table_full_3_test.csv")
write.csv(merged, 
          "~/Desktop/rodrigo/public-clonotype/Analysis/2020-08-11_public_table_summary.csv",
          row.names = F)



#redoing tree and trying to select the clones found in more than one animal
x_m_summary_select <- x_m_summary %>% 
  group_by(clonoquery, name, CDR3_aa, VDJ_nt) %>%
  summarise()

public <- x_m_basic_table %>% filter(nchar(animals) > 3)

#write the code for the tree for clones expanded
x_m_summary_select








