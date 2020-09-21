library(igraph)
library(data.table)
library(ape)
library(Biostrings)
library(msa)

#parallel processing
require(doParallel)
cores <- 2
system <- Sys.info()['sysname']
cl <- NULL
if (system == 'Windows') {
  cl <- makeCluster(getOption('cl.cores', cores))
  registerDoParallel(cl)
  registerDoSEQ()
  on.exit(stopCluster(cl))
} else {
  options('mc.cores' = cores)
  registerDoParallel(cores)
}


#igraph test
ls <- list.files("~/Desktop/rodrigo/igraph/Clonoquery full rep/", full.names = T)
cl_f <- lapply(ls,fread,fill = T, blank.lines.skip = T, 
            stringsAsFactors = F,
            header = T)
#remove empty rows with 0 clones 
cl_omit <- lapply(cl_f, na.omit)


#write.fast
x <- readDNAStringSet("~/Desktop/rodrigo/igraph/E11_t.fasta",
                      format="fasta")
names(x) <- cl_omit[[1]][["name"]]
#ignore-testing
paste(cl_omit[[1]][["name"]], str_count(x),sep = "_")

if (duplicate(x)){
  names(x) <- paste0(x,"_",)
}
  (duplicated(x), paste, paste0(x,)
make.unique(names(x))
 duplicated(x)      
       unique(x)

x_edited <- ifelse(x)

#write fasta
writeXStringSet(x,
            file = "~/Desktop/rodrigo/igraph/E11_t.fasta")

#align

aligned <- msaMuscle(x)


#distance
dist <- readDNAMultipleAlignment("~/Desktop/rodrigo/igraph/e11_aligned.fasta")
dist_b <- as.DNAbin(dist)
dist_c <- dist.dna(dist_b)

str(dist_c)
x <- matrix(dist_c, nrow = 3,ncol=3)
View(x)

x
g <- graph.adjacency(
  as.matrix(dist_c),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

corx <- as.dist(cor(t(dist_c)))
g <- graph.adjacency(
  as.matrix(as.dist(cor(t(x), method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)
g1 <- plot(g)
g1


