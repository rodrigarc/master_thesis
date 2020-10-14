

sc_df <- read.table("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/sc-v4_filtered-files/2020-08-10_sc-combined-filtered.tab",
           header = T, sep = "\t", quote = "")
sc_df <- sc_df %>% mutate(sequence_1 = gsub("_HC_PostF|_HC_PreF+|_HC_PreF|_HC_DP|_HC_DP1|_HC_DP2|+",
                             "",sc_df$name)) 
df_comb <- df_comb %>% mutate(sequence_1 = gsub("_HC_PostF|_HC_PreF+|_HC_PreF|_HC_DP|_HC_DP1|_HC_DP2|+",
                                              "", df_comb$sequence)) 

ptn <- paste(sc_df$sequence_1,collapse="\\b|\\b")
grep(paste0("\\b", ptn,"\\b"), 
     df_comb$sequence_1, value=TRUE)

matches <- unique(grep(ptn, 
                       df_comb$sequence_1, value=TRUE))

df_comb_cdrs <- df_comb[df_comb$sequence_1 %in% matches]

ptn
