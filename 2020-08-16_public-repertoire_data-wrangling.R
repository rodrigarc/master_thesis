library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(stringr)

getwd()
setwd("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis")

ls <- list.files("v3", full.names = T)
ls <- grep("2020-08-13" ,ls,value = T)
ls <- grep("csv" ,ls,value = T)
ls <- lapply(ls,read.csv)
names(ls) <- c("B1","B2","IgM")
df <- rbindlist(ls,idcol=T, fill = T)


#new table
ndf <- select(df, sequence, clone,.id)
ndf_e <- ndf %>% mutate(IgM = ifelse(grepl("IgM", .id), 1, 0),
                        B1 = ifelse(grepl("B1", .id), 1, 0),
                        B2 = ifelse(grepl("B2", .id), 1, 0)) %>%
  select(-.id) %>% group_by(clone, sequence) %>% summarise(IgM = sum(IgM),
                                                           B1 = sum(B1),
                                                           B2 = sum(B2))
ndf_s <- ndf_e %>% separate_rows(sequence, sep = ", ") %>%
  rename(well = sequence)
ndf_s$well <- gsub("+","",ndf_s$well, fixed = T)

#write.csv(ndf_s, "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v3/2020-08-16_public-table-combined.csv",
row.names = F, quote = F)

#combining animals
rem_dup.one <- function(x){
  paste(unique(toupper(trimws(unlist(strsplit(x,split="(?!')[ [:punct:]]",fixed=F,perl=T))))),collapse = ", ")
}
rem_dup.vector <- Vectorize(rem_dup.one,USE.NAMES = F)

df_f <- df %>% 
  mutate(combined_animals = paste0(animals_rep_seq, sep = ", ", animals_sc)) 

df_f$combined_animals <- rem_dup.vector(df_f$combined_animals)
df_f$combined_animals <- gsub(", , ", ", ",df_f$combined_animals)
df_f$.id <- gsub("IgM", "PV", df_f$.id)
df_f$sequence <- sub("+", "",df_f$sequence)
df_comb <- df_f %>% mutate(number_animals_combined = str_count(combined_animals, "E"))

#groups based on NP and SOL
groups <- data.frame(NP = c("E11", "E16", "E17", "E23", "E24"),
                     SOL = c("E12", "E14", "E18", "E21", NA))
ff = function(x, patterns, replacements = patterns, fill = NA, ...)
{
  stopifnot(length(patterns) == length(replacements))
  
  ans = rep_len(as.character(fill), length(x))    
  empty = seq_along(x)
  
  for(i in seq_along(patterns)) {
    greps = grepl(patterns[[i]], x[empty], ...)
    ans[empty[greps]] = replacements[[i]]  
    empty = empty[!greps]
  }
  
  return(ans)
}
df_comb$NP <- ff(df_comb$combined_animals, 
                     patterns = c("E11", "E16", "E17", "E23", "E24"), 
                     replacements = c("NP","NP","NP","NP","NP"))
df_comb$SOL <- ff(df_comb$combined_animals, 
                  patterns = c("E12", "E14", "E18", "E21"), 
                  replacements = c("SOL","SOL","SOL","SOL"))

df_comb <- df_comb %>%
  mutate(Groups = ifelse(!is.na(df_comb$NP) & !is.na(df_comb$SOL),
                         "Both",
                         ifelse(!is.na(df_comb$NP) & is.na(df_comb$SOL),
                                "NP",
                                ifelse(is.na(df_comb$NP) & !is.na(df_comb$SOL),
                                       "SOL", ""))))
#dataframe with combined animals
new_df <- df_f %>% group_by(clone) %>% 
  summarise(combined = paste(combined_animals, collapse = ", ")) %>%
  group_by(clone, combined) 
new_df$combined <- rem_dup.vector(new_df$combined)
new_df$combined <- gsub(", , ", ", ",new_df$combined)
new_df$animals_sc_rep_seq <- str_count(new_df$combined, "E")
new_df <- new_df %>% arrange(desc(animals_sc_rep_seq))

#write.csv(new_df,
          "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v3/2020-08-18_public-clones-sc-rep-sc.csv",
          row.names = F)


#testing plots
set.seed(10)
df_comb$.id <- factor(df_comb$.id ,levels = c("PV", "B1", "B2"))
g1 <- df_comb %>%
  ggplot(aes(x = .id, y = number_animals_combined, size = size_rep_seq,
             color = Groups)) +
  geom_point(position = "jitter") +
  scale_color_manual(values = c("#107F01","#065386", "#5C003C"))

g2 <- g1 + theme_cowplot() + xlab("") +ylab("Number of animals") +
  labs(size = "Clone size") + theme(legend.position = "right") +
  scale_size_continuous(breaks = c (10, 100,1000,10000), limits = c(0,10000))+
  scale_y_continuous(breaks = c(2,3,4,5,6,7,8,9)) 

g3 <- g2 + geom_vline(xintercept=c(-3,1.5), linetype="dotted") +
  geom_vline(xintercept=c(-3,2.5), linetype="dotted")

g3


#ggsave("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v3/2020-08-17_public-clones-IgM-B1-B2.tiff",
       height = 3.5,
       width = 5,
       dpi = 200,
       bg = "transparent")

minors<-seq(1.5,71,by=1)
minors
df_comb %>% ggplot(aes(x = clone, 
                       size = size_rep_seq, 
                       y = number_animals_combined,
                       fill = .id)) +
  geom_point(position = position_dodge(width = 0.5), 
             shape = 21, color = "black",stroke = 0.5) +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90,vjust = .5, hjust=1),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank())+
  scale_color_viridis_d()+
  labs(
    x = "Clones",
    y = "Number of animals",
    fill = "Timepoints",
    size = "Clone size"
  )+
  geom_vline(mapping=NULL, xintercept = minors,
             colour='grey50', linetype = "dotted")+
  scale_size_continuous(breaks = c (10, 100,1000,10000), limits = c(0,10000))

#ggsave("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v3/2020-08-17_all-public-clones-IgM-B1-B2.tiff",
       height = 7,
       width = 12,
       dpi = 150,
       bg = "transparent")
