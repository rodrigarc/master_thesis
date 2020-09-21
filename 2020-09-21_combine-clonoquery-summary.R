library(data.table)
library(dplyr)

###create function to read and combine data 
##data = clonoquery_summary located in a folder containing all the txt files
##all the animals files should be in the folowing order inside the folder:
##c(E11", "E12", "E14", "E16", "E17", "E18", "E21", "E23", "E24" 

#x = file location with directory containing summary clonoqueries
.clono_combine <- function(x){
  ls <- list.files(x, full.names = T)
  print(ls)
  names(ls) <- c("E11", "E12", "E14", "E16", "E17", "E18", "E21", "E23", "E24")
  clonoquery_summary <- lapply(ls, fread, fill=TRUE)
  clonoquery_summary <- rbindlist(clonoquery_summary, idcol = TRUE)
  clonoquery_summary <- clonoquery_summary[,-c(18,19)] %>% rename(ID = .id)
  return(clonoquery_summary)
}

###### COMBINING CLONOQUERY SUMMARY v4_NGS TR1 #######
clonoquery_summary_B2 <- .clono_combine(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_NGS TR1/Clonoquery summary")

###### COMBINING CLONOQUERY SUMMARY v4_B1-D14 #######
clonoquery_summary_B1 <- .clono_combine(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_B1-D14/clonoquery summary")

###### COMBINING CLONOQUERY SUMMARY igM #######
clonoquery_summary_igm <- .clono_combine(
  "~/Box Sync/RSV NGS/v4_Analysis/v4_igm_clonoquery/clonoquery summary")


combined <- list(clonoquery_summary_igm, clonoquery_summary_B1, clonoquery_summary_B2)
names(combined) <- c("PV", "B1", "B2") 
combined <- rbindlist(combined, idcol = T)  
combined <- combined %>% rename(time_point = .id)

write.csv(combined, "~/Desktop/RSV/2020-09-21_clonoquery_summary.csv",row.names = F)
  
  
  
  