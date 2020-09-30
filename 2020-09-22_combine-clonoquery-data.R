library(data.table)
library(dplyr)

###### Create Function to read and combine data from clonoquery summary
##data = clonoquery_summary located in a folder containing all the txt files per animal
##all the animals files should be in the folowing order inside the folder:
##c(E11", "E12", "E14", "E16", "E17", "E18", "E21", "E23", "E24" 
##both functions print the list of files inside the directory you selected as input
##always check the list to see if it is in the right order of animals and correct files

#x = file location with directory containing summary clonoqueries
.clono_summary <- function(x){
  ls <- list.files(x, full.names = T)
  if(length(ls) != 9) {
    stop("Folder not containing 9 files, please check the number of files or change the function to allow a different number of files.")
  }
  if(any(!grepl("summary",ls))){
    stop("Files in the folder are not identified as a summary files.")
  }
  else{
    print(ls)
    names(ls) <- c("E11", "E12", "E14", "E16", "E17", "E18", "E21", "E23", "E24")
    clonoquery_summary <- lapply(ls, fread, fill=TRUE)
    clonoquery_summary <- rbindlist(clonoquery_summary, idcol = TRUE)
    clonoquery_summary <- clonoquery_summary[,-c(18,19)] %>% rename(ID = .id)
  }
  return(clonoquery_summary)
}

###### COMBINING CLONOQUERY SUMMARY v4_NGS TR1 #######
clonoquery_summary_B2 <- .clono_summary(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clonoquery summary")

###### COMBINING CLONOQUERY SUMMARY v4_B1-D14 #######
clonoquery_summary_B1 <- .clono_summary(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_B1-D14/clonoquery summary")

###### COMBINING CLONOQUERY SUMMARY igM #######
clonoquery_summary_igm <- .clono_summary(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_igm_clonoquery/clonoquery summary")

##### Combining all the clonoquery summary files into a single file #######
## output as CSV file with a new column for time_point

combined_summary <- list(clonoquery_summary_igm, clonoquery_summary_B1, clonoquery_summary_B2)
names(combined_summary) <- c("PV", "B1", "B2") 
combined_summary  <- rbindlist(combined_summary, idcol = T)  
combined_summary  <- combined_summary %>% rename(time_point = .id)

write.csv(combined_summary, "~/Desktop/RSV/2020-09-21_clonoquery_summary.csv",row.names = F)
  
###### Create Function to read and combine data from clonoquery full
##argument is the same as in 
.clono_full <- function(z){
ls <- list.files(z, full.names = T)
  if(length(ls) != 9) {
    stop("Folder not containing 9 files, please check the number of files or change the function to allow a different number of files.")
  }
  if(any(!grepl("full",ls))){
    stop("Files in the folder are not identified as a full file in its name.")
  }
  else{
    print(ls)
    x <- lapply(ls, fread, fill = T)
    names(x) <- c("E11", "E12", "E14", "E16", "E17", "E18", "E21", "E23", "E24")
    x <- lapply(x, function (x) x[!apply(is.na(x) | x == "", 1, all),])       
    x_all <- do.call(rbind, x)
#getting the clonoquery named properly, all clones queries, aka full table
    x_m_full <- x_all %>% 
      mutate(grp = cumsum(grepl("Query",name))) %>%
      group_by(grp) %>% 
      mutate(clonoquery = ifelse(row_number() > 1, first(name), name)) %>%
      ungroup() %>% 
      select(-c(grp)) %>%
      relocate(clonoquery)
    }
return(x_m_full)
}

###### COMBINING CLONOQUERY FULL v4_NGS TR1 #######
clonoquery_full_B2 <- .clono_full(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clonoquery full")

###### COMBINING CLONOQUERY FULL v4_B1-D14 #######
clonoquery_full_B1 <- .clono_full(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_B1-D14/clonoquery full")

###### COMBINING CLONOQUERY FULL igM #######
clonoquery_full_igm <- .clono_full(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_igm_clonoquery/clonoquery full/")

##### Combining all the clonoquery full files into a single file #######
## output as CSV file with a new column for time_point

combined_full <- list(clonoquery_full_igm, clonoquery_full_B1, clonoquery_full_B2)
names(combined_full) <- c("PV", "B1", "B2") 
combined_full  <- rbindlist(combined_full, idcol = T, fill = T)  
combined_full  <- combined_full %>% rename(time_point = .id)

write.csv(combined_full, "~/Desktop/RSV/2020-09-21_clonoquery_full.csv",row.names = F)
