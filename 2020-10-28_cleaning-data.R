######## Processing small data for Ben Murrel analysis ########
library(dplyr)
library(data.table)

#read filtered tab file containing the igdiscover output from the bulk B2 of animal E16 with bulk germline
df <- fread("~/Box Sync/Ben Murrel presentation/E16_B2_filtered.tab.gz", quote = F)

#save a smaller version of the dataframe compressed in case you need it
df_sel <- df %>% select(name,V_gene, J_gene, VDJ_nt, VDJ_aa) 
fwrite(x = df_sel, file = "~/Box Sync/Ben Murrel presentation/E16_B2_clean.tab.gz",quote = F )

#save tab for nt and aa fasta file conversion
df_sel_nt <- df %>% select(name, VDJ_nt) 
fwrite(x = df_sel, file = "~/Box Sync/Ben Murrel presentation/E16_B2_clean_nt.tab.gz",quote = F )

df_sel_aa <- df %>% select(name, VDJ_aa) 
fwrite(x = df_sel, file = "~/Box Sync/Ben Murrel presentation/E16_B2_clean_aa.tab.gz",quote = F )


####### DO NOT USE IT, USES TOO MUCH MEMORY ###save it as fasta as well
df_fasta <- df %>% select(name, VDJ_nt) 
dataframe2fas(df_fasta, "~/Box Sync/Ben Murrel presentation/E16_B2_clean.fasta")
