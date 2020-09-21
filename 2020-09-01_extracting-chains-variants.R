library(Biostrings)
library(openxlsx)

#did not use it, because I have done it before, however the code is good
#to combine many fasta files
ls <- list.files("~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined", 
                 full.names = T)[-10]
db <- lapply(ls, readDNAStringSet)
db_c <- do.call(c, db)

#using files
HJ <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/Summary germlines/HJ aa.fasta")
HV <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/combined/v_aa.fasta")
LJ <-  readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/Summary germlines/LJ aa.fasta")
LV <-  readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_Germlines/Summary germlines/LV aa.fasta")

#extracting first 8 AA from J gene and last 8 AA from V gene
HV_start <- BStringSet(HV, start=1, end = 8)
HJ_end <- BStringSet(HJ, start=-8, end = -1)
LV_start <- BStringSet(LV, start=1, end = 8)
LJ_end <- BStringSet(LJ, start=-8, end = -1)
#selecting unique variants
HV_unique <- unique(HV_start)
HJ_unique <- unique(HJ_end)
LV_unique <- unique(LV_start)
LJ_unique <- unique(LJ_end)

writeXStringSet(HJ_unique, 
                file = "~/Box Sync/RSV NGS/mAb_QC/HJ_unique.fasta")
writeXStringSet(HV_unique, 
                file = "~/Box Sync/RSV NGS/mAb_QC/HV_unique.fasta")
writeXStringSet(LV_unique, 
                file = "~/Box Sync/RSV NGS/mAb_QC/LV_unique.fasta")
writeXStringSet(LJ_unique, 
                file = "~/Box Sync/RSV NGS/mAb_QC/LJ_unique.fasta")

#matching patterns
HC_all <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/AA sequences/LOR-HC.fasta")
LC_all <- readAAStringSet(
  "~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/AA sequences/LOR-LC.fasta")

BStringSet(HC_all, start=1, end = 8) %in% HV_unique 
BStringSet(HC_all, start=-8, end = -1) %in% HJ_unique
BStringSet(LC_all, start=1, end = 8) %in% LV_unique 
BStringSet(LC_all, start=-8, end = -1) %in% LJ_unique



#testing as functions
#checking if unique patterns are not found and flag as 1
HC_check <- function(HC_seqs, HV_pattern, HJ_pattern){
HC<- data.frame(name = names(HC_seqs))
HC$mismatch_start_HC <- ifelse(
  !BStringSet(HC_seqs, start=1, end = 8) %in% HV_pattern,
         1, 0)
HC$mismatch_end_HC <- ifelse(
  !BStringSet(HC_seqs, start=-8, end = -1) %in% HJ_pattern,
  1, 0)
return(HC)
print(HC)
}


HC_check(HC_all, HV_unique, HJ_unique)

#test function for evalue and coverage
function(x){
  ifelse (x$V_coverage_HC < 99, 1 , 0)
  ifelse (x$V_coverage_LC < 99, 1 , 0)
  ifelse (x$V_evalue_HC < 2.52e-114, 1 , 0)
  ifelse (x$V_evalue_LC < 3.51e-113, 1 , 0)
  ifelse (x$J_coverage_HC < 85.4, 1 , 0)
  ifelse (x$J_coverage_LC < 92.1, 1 , 0)
  ifelse (x$J_evalue_HC < 4.49e-17, 1 , 0)
  ifelse (x$J_evalue_LC < 9.61e-14, 1 , 0)
}
df$mAb.ID == ab$name

#combined function
df <- read.xlsx(
  "~/Box Sync/RSV NGS/v4_Analysis/LOR mAbs/200810 mAb clones.xlsx",
  sheet = 2)

mab_check <- function(LC_seqs, LV_pattern, LJ_pattern,
                      HC_seqs, HV_pattern, HJ_pattern, df_scores){
  if (missing(LC_seqs)) {
    warning("Please add all the 7 necessary arguments")
  }
  if (any(names(HC_seqs) != names(LC_seqs))) {
    warning("Heavy chain and light chain names do not match.")
  }
    mab<- data.frame(name = names(HC_seqs))
    mab$mismatch_start_LC <- ifelse(
      !BStringSet(LC_seqs, start=1, end = 8) %in% LV_pattern,
      1, 0)
    mab$mismatch_end_LC <- ifelse(
      !BStringSet(LC_seqs, start=-8, end = -1) %in% LJ_pattern,
      1, 0)
    mab$mismatch_start_HC <- ifelse(
      !BStringSet(HC_seqs, start=1, end = 8) %in% HV_pattern,
      1, 0)
    mab$mismatch_end_HC <- ifelse(
      !BStringSet(HC_seqs, start=-8, end = -1) %in% HJ_pattern,
      1, 0)
    mab$HV_coverage <- ifelse(df_scores$HV.coverage < 99, 1 , 0)
    mab$HV_evalue <- ifelse (df_scores$HV.evalue < 2.52e-114, 1 , 0)
    mab$HJ_coverage <- ifelse (df_scores$HJ.coverage < 85.4, 1 , 0)
    mab$HJ_evalue <- ifelse (df_scores$HJ.evalue < 4.49e-17, 1 , 0)
    mab$LV_coverage <- ifelse (df_scores$LV.coverage < 99, 1 , 0)
    mab$LV_evalue <- ifelse (df_scores$LV.evalue < 3.51e-113, 1 , 0)
    mab$LJ_coverage <- ifelse (df_scores$LJ.coverage < 92.1, 1 , 0)
    mab$LJ_evalue <- ifelse (df_scores$LJ.evalue < 9.61e-14, 1 , 0)
    mab <- cbind(mab, total = rowSums(mab[,2:ncol(mab)]))
  return(mab)
  print(mab)
}

ab <- mab_check(LC_all, LV_unique, LJ_unique,
          HC_all, HV_unique, HJ_unique, df)

kableExtra::kable(ab)
