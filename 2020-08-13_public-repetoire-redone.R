#rename single-cell clontypes by number
#query using clonotypes file (which is filtered tab)



#rethink how to calculate the distance between each clonotype and 
#download clonoquery TR1_v4 for every animal and then calculate
#between every animal
#plot using igraph

library(data.table)
library(dplyr)

x <- fread("~/Desktop/rodrigo/public-clonotype/Public_clonotypes/Public_clonotypes_full.txt",
           fill = T)
x <- x[!apply(is.na(x) | x == "", 1, all),]


x$name <- make.unique(x$name)

fwrite(x,"~/Desktop/rodrigo/public-clonotype/Public_clonotypes/Public_clonotypes_edited_full.txt",
       sep = "\t", quote = F)

#after rerunning clonoquery with a unique name to each clone
ls <- list.files("~/Desktop/rodrigo/public-clonotype/Public_clonoquery_3/B1", 
                full.names = T,
                pattern = "summary")
x <- lapply(ls, fread, fill = T)
names(x) <- ls
x <- lapply(x, function (x) x[!apply(is.na(x) | x == "", 1, all),])       
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
knitr::kable(x_m_basic_table) %>%
kable_styling(bootstrap_options = "striped", full_width = F)
  


#csv files from dfs
write.csv(x_m_full, 
          "~/Desktop/rodrigo/public-clonotype/Analysis/Public_table_full.csv")
write.csv(x_m_summary, 
          "~/Desktop/rodrigo/public-clonotype/Analysis/Public_table_summary.csv")

#redoing tree and trying to select the clones found in more than one animal
x_m_summary_select <- x_m_summary %>% 
  group_by(clonoquery, name, CDR3_aa, VDJ_nt) %>%
  summarise()

public <- x_m_basic_table %>% filter(nchar(animals) > 3)

#write the code for the tree for clones expanded
x_m_summary_select








