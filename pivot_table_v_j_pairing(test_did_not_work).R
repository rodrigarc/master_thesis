library(openxlsx)
library(dplyr)

ls <- read.xlsx("~/Desktop/rodrigo/V-J Pairing /List of V and J genes.xlsx")
ls <- ls %>% dplyr::rename(V_gene = Unique.Genes, J_gene = J.Genes)
vj <- read.xlsx("~/Desktop/rodrigo/V-J Pairing /VJ frequencies_TR1.xlsx")
vj <- vj %>% dplyr::rename(J_gene = J.gene)

library(dplyr)
library(tidyr)

ls <- data.frame(A = c("ABC", "DEF", "GHI","XYZ"),
                 B = c("KLM","MNO", "",""))

df <- data.frame(ID = c(1,1,2,2),
                 S = c("x","y","y","z"),
                 A = c("ABC","DEF","DEF","XYZ"), 
                 B = c("KLM","MNO","XYZ","MNO"), 
                 C = c(100,150,2,10))

36-24
df %>%
  complete(S, A = unique(ls$A), B, fill = list(C = 0)) %>%
  group_by(S) %>%
  fill(ID, .direction = "downup") %>% View()

df %>% group_by(ID, S) %>% 
  ifelse(ls$A %in% df$A & ls$B %in% df$B, "",add_row(ID = df$ID,
                                                     S = df$S,
                                                     A = ls$A,
                                                     B = ls$B,
                                                     C = 0))
lsb <- subset(ls,ls$V_gene != "")
ls2 <- expand.grid(Specificity=vj$Specificity, V_gene=ls$V_gene, J_gene=lsb$J_gene)
ls3 <- expand.grid(ID=vj$ID, V_gene=ls$V_gene, J_gene=lsb$J_gene)
ls4 <- cbind(ID=ls3$ID,ls2)

lsa <- mutate(ls4, Size =ifelse((ls4$ID==vj$ID & ls4$Specificity==vj$Specificity & 
                                   ls4$V_gene==vj$V_gene & ls4$J_gene==vj$J_gene) , vj$Size, 0))
View(lsa)

test <- vj %>%
  complete(Specificity, 
           V_gene = unique(ls$V_gene), J_gene, 
           fill = list(Size = 0)) %>%
  group_by(Specificity) %>%
  fill(ID, .direction = "downup") %>%
  distinct(.keep_all = TRUE) %>%
  select(ID, Specificity, V_gene, J_gene, Size)

vj_sum_test <- test  %>% 
  filter(Specificity != "PostF") %>%
  group_by(ID, Specificity, V_gene,J_gene) %>% 
  mutate(n = sum(Size)) %>%
  ungroup() %>% group_by(ID, Specificity) %>% 
  group_split() 

vj_sum_test_2 <- lapply(vj_sum_test, 
                        function (x) 
                          complete(x,
                                   Specificity, 
                                   V_gene = unique(ls$V_gene), 
                                   J_gene = unique(na.omit(ls$J_gene)), 
                                   fill = list(Size = 0, n = 0)) %>%
                          group_by(Specificity) %>%
                          fill(ID, .direction = "downup") %>%
                          distinct(.keep_all = TRUE) %>% 
                          group_by(ID, Specificity, V_gene,J_gene) %>% 
                          mutate(n = sum(Size)) %>%
                          select(ID, Specificity, V_gene, J_gene,n))

vj_total <- vj %>% group_by(ID,V_gene,J.gene) %>% 
  mutate(n = sum(Size)) %>%
  ungroup() %>% group_by(ID) %>% 
  group_split() 


?grepl
