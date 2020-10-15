library(data.table)
library(seqRFLP)
library(Biostrings)
library(dplyr)
library(stringr)

##creating the df_comb and combined full (combining que full clonoquery)
setwd("~/Box Sync/public_repertoire/results/2020-09-13/")
ls <- list.files("public-tables", full.names = T)
ls <- grep("2020-09-13" ,ls,value = T)
ls <- grep("summary" ,ls,value = T)
ls
ls <- lapply(ls,read.csv)
names(ls) <- c("B1","B2","IgM")
df <- rbindlist(ls,idcol=T, fill = T)

lc_all <- openxlsx::read.xlsx("~/Box Sync/RSV NGS/v4_Analysis/LCs/200828_LC_public.xlsx", sheet = 2)


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

df_comb_jit <- transform(df_comb, dmeasure =  jitter(as.numeric(as.factor(.id)), 2),
                         jvalue = jitter(number_animals_combined, amount = .2) )
df_comb_jit$sequence <- gsub("_HC_PostF|_HC_PreF+|_HC_PreF|_HC_DP|_HC_DP1|_HC_DP2|+",
                             "",df_comb_jit$sequence)
ptn <- paste(lc_all$well,collapse="\\b|\\b")
matches <- unique(grep(paste0("\\b", ptn,"\\b"), 
                       df_comb_jit$sequence, value=TRUE))
df_comb_jit_with_LC <- df_comb_jit[df_comb_jit$sequence %in% matches]
##### combinin clonoquery function
.clono_full <- function(z){
  ls <- list.files(z, full.names = T)
  ls <- grep("full.txt", ls, value = T)
  if(length(ls) != 9) {
    stop("Folder not containing 9 files, please check the number of files or change the function to allow a different number of files.")
  }
  if(any(!grepl("full",ls))){
    stop("Files in the folder are not identified as a full file in its name.")
  }
  else{
    x <- lapply(ls, fread, fill = T)
    names(x) <- c("E11", "E12", "E14", "E16", "E17", "E18", "E21", "E23", "E24")
    x <- lapply(x, function (x) x[!apply(is.na(x) | x == "", 1, all),])       
    x_all <- rbindlist(x, idcol = T)
    x_all <- x_all %>% rename(ID = .id)
    #getting the clonoquery named properly, all clones queries, aka full table
    x_m_full <- x_all %>% 
      mutate(grp = cumsum(grepl("Query",name))) %>%
      group_by(grp) %>% 
      mutate(clonoquery = ifelse(row_number() > 1, first(name), name)) %>%
      ungroup() %>% 
      select(-c(grp)) %>%
      relocate(clonoquery) %>%
      mutate(clone = substr(clonoquery, 10,14))
  }
  return(x_m_full)
}

###### COMBINING CLONOQUERY FULL v4_NGS TR1 #######
clonoquery_full_B2 <- .clono_full(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/clonoquery/B2")

###### COMBINING CLONOQUERY FULL v4_B1-D14 #######
clonoquery_full_B1 <- .clono_full(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/clonoquery/B1")

###### COMBINING CLONOQUERY FULL igM #######
clonoquery_full_igm <- .clono_full(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-09-13/clonoquery/IgM/")

##### Combining all the clonoquery full files into a single file #######
## output as CSV file with a new column for time_point

combined_full <- list(clonoquery_full_igm, clonoquery_full_B1, clonoquery_full_B2)
names(combined_full) <- c("PV", "B1", "B2") 
combined_full  <- rbindlist(combined_full, idcol = T, fill = T)  
combined_full  <- combined_full %>% 
  rename(time_point = .id) %>% 
  filter(count > 0) %>%
  mutate(full_ID = paste0(clone, ID, time_point,))
  


######selecting public clones #######
df_comb_selection <- df_comb_jit_with_LC %>% filter(number_animals_combined >= 5)
p_clone_selection <- unique(df_comb_selection$clone)
p_clone_sequences <- combined_full[combined_full$clone %in% p_clone_selection]

##### reading files and deleting * character ######
p_nt <- p_clone_sequences %>% select(name, VDJ_nt, clone) 
p_nt$VDJ_nt <- gsub("\\*","", p_nt$VDJ_nt)
w <- "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/not_dedup/"

p_aa <- p_clone_sequences %>% select(name, VDJ_aa, clone)
p_aa$VDJ_aa <- gsub("\\*","", p_aa$VDJ_aa)

#### loops for creating the fasta files out of the dataframes per clone"
for (i in unique(p_nt$clone)){
dataframe2fas(p_nt[p_nt$clone == i, c("name","VDJ_nt")], paste0( w, i,"_public_clones_5animals_nt.fasta"))
}

for (i in unique(p_aa$clone)){
  dataframe2fas(p_aa[p_aa$clone == i, c("name","VDJ_aa")], paste0( w, i,"_public_clones_5animals_aa.fasta"))
}

#### you can check the tables from this to check the number of repeated seqs

ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/not_dedup/",full.names = T)
ls_aa <- grep("aa", ls, value = T)
ls_nt <- grep("nt", ls, value = T)


#### reading fasta files and deduplicating with loops
fasta_aa <- sapply(ls_aa,readBStringSet,simplify = FALSE,USE.NAMES = TRUE)
fasta_nt <- sapply(ls_nt,readBStringSet,simplify = FALSE,USE.NAMES = TRUE)



fasta_aa_u <- lapply(fasta_aa, unique)
names(fasta_aa_u) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")
fasta_nt_u <- lapply(fasta_nt, unique)
names(fasta_nt_u) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")

fasta_nt_d <- lapply(fasta_nt, function(x) x[duplicated(x)])
names(fasta_nt_d) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")

duplicates <- list()
counts_l <- list()
for (i in names(fasta_nt_d)){
  duplicates <- match(fasta_nt_d[[i]], fasta_nt_u[[i]])
  counts_l[[i]] <- data.frame(table(duplicates))
  n <- as.numeric(counts_l[[i]]$duplicates)
  counts_l[[i]][["names"]] <- names(fasta_nt_u[[i]][n])
}

z <- "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/dedup/"

for (i in names(fasta_aa_u)){
  writeXStringSet(fasta_aa_u[[i]], paste0(z, i, "_dedup_public_clones_5animals_aa.fasta"))
}

for (i in names(fasta_nt_u)){
  writeXStringSet(fasta_nt_u[[i]], paste0(z, i, "_dedup_public_clones_5animals_nt.fasta"))
}

####### repeating the whole process for NANOPARTICLE  GROUP (NP) only ##### 

df_comb_selection_np <- df_comb_jit_with_LC %>% filter(Groups == "NP")
p_clone_selection_np <- unique(df_comb_selection_np$clone)
p_clone_sequences_np <- combined_full[combined_full$clone %in% p_clone_selection_np]

##### reading files and deleting * character ######
p_nt_np <- p_clone_sequences_np %>% select(name, VDJ_nt, clone) 
p_nt_np$VDJ_nt <- gsub("\\*","", p_nt_np$VDJ_nt)
w <- "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/NP/"

p_aa_np <- p_clone_sequences_np %>% select(name, VDJ_aa, clone)
p_aa_np$VDJ_aa <- gsub("\\*","", p_aa_np$VDJ_aa)

#### loops for creating the fasta files out of the dataframes per clone"
for (i in unique(p_nt_np$clone)){
  dataframe2fas(p_nt_np[p_nt_np$clone == i, c("name","VDJ_nt")], paste0( w, i,"_public_clones_np_nt.fasta"))
}

for (i in unique(p_aa_np$clone)){
  dataframe2fas(p_aa_np[p_aa_np$clone == i, c("name","VDJ_aa")], paste0( w, i,"_public_clones_np_aa.fasta"))
}

#### you can check the tables from this to check the number of repeated seqs

ls_np <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/NP/",full.names = T)
ls_aa_np <- grep("aa", ls, value = T)
ls_nt_np <- grep("nt", ls, value = T)


#### reading fasta files and deduplicating with loops
fasta_aa_np <- sapply(ls_aa_np,readBStringSet,simplify = FALSE,USE.NAMES = TRUE)
fasta_nt_np <- sapply(ls_nt_np,readBStringSet,simplify = FALSE,USE.NAMES = TRUE)

fasta_aa_u_np <- lapply(fasta_aa_np, unique)
names(fasta_aa_u_np) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")
fasta_nt_u_np <- lapply(fasta_nt_u_np, unique)
names(fasta_nt_u) <- c("P0049", "P0108", "P0254", "P0733", "P0741", "P0867", "P0913", "P1039", "P1123", "P2029")


z <- "~/Box Sync/RSV NGS/v4_Analysis/v4_public-repertoire/analysis/v4/results/2020-10-13/NP/"

for (i in names(fasta_aa_u)){
  writeXStringSet(fasta_aa_u[[i]], paste0(z, i, "_dedup_public_clones_NP_aa.fasta"))
}

for (i in names(fasta_nt_u)){
  writeXStringSet(fasta_nt_u[[i]], paste0(z, i, "_dedup_public_clones_NP_nt.fasta"))
}

