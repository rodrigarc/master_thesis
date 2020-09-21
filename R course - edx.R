install.packages("dslabs")
library(dslabs)
data(heights)
str(heights)

heights[777,]
heights[777,1]
max(heights$height)
which.min(heights$height)

mean(heights$height)

nrow(filter(heights, sex == "Male"))
nrow(heights)
812/1050

heights %>% filter(sex == "Female" & height >= 78 ) 

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install(c("genefilter","geneplotter"))
library(BSgenome.Hsapiens.UCSC.hg19)
library(genefilter)
library(geneplotter)
yes

BiocManager::install(c("genefu",
                       "COPDSexualDimorphism.data",
                       "GenomicRanges",
                       "gwascat",
                       "hgu133a.db",
                       "genomicsclass/tissuesGeneExpression"))
BiocManager::install("GenomeInfoDbData")
BiocManager::install("rtracklayer")
BiocManager::install("SummarizedExperiment")

library("rtracklayer")
library(genefu)
data(sig.gene70)
dim(sig.gene70)
head(sig.gene70)[,1:6]
sum(is.na(sig.gene70$NCBI.gene.symbol))
sig.gene70 %>% filter(Description == "cyclin E2")

grep("kinase", sig.gene70$Description)

library(COPDSexualDimorphism.data)
data(lgrc.expr.meta)
table(expr.meta$GENDER)
summary(expr.meta$pkyrs)
qqnorm(expr.meta$pkyrs)
boxplot(pkyrs~gender, data=expr.meta)


library(GenomicRanges)
BiocManager::install("gwascat")
library(gwascat)
data(ebicat37)
ebicat37
sort(table(ebicat37$CHR_ID),decreasing=TRUE)

x <- sort(table(mcols(ebicat37)[,"DISEASE/TRAIT"]), decreasing=TRUE)
x[1]

#day 2020-07-23
#Bioconductor and GitHub packages
BiocManager::install(c("Homo.sapiens",
                       "GenomicFeatures",
                       "genomicsclass/ERBS",
                       "genomicsclass/ph525x"))

library(ERBS)
data(HepG2)
class(HepG2)
HepG2

data("GM12878")
seqnames(HepG2)
chr <- seqnames(HepG2)
HepG2[chr=="chr20",]
table(chr)[1:24] 

x <- mcols(HepG2)
median(x$signalValue)
x
which.max(x$signalValue)

HepG2[120]
hist(table(chr))
median(table(chr))


median(width(HepG2))
hist(width(HepG2),nclass=25)

library(IRanges)
x <- IRanges(101,200)
x*2
narrow(x, start=20)
x+25

x <- IRanges(start=c(1,11,21), end=c(3,15,27))
x
sum(width(x))

x <- IRanges(start=c(101,106,201,211,221,301,306,311,351,361,401,411,501), 
             end=c(150,160,210,270,225,310,310,330,390,380,415,470,510))

library(ph525x)
plotRanges(x)
sum(width(gaps(x)))

disjoin(x)

par(mfrow=c(2,1))
plotRanges(x,xlim=c(0,600))
plotRanges(resize(x,50), xlim=c(0,600))
?resize
resize(x,1)
IRanges(x)

library(IRanges)
library(GenomicRanges)
library(ph525x)

x = GRanges("chr1", IRanges(c(1,101),c(50,150)), strand=c("+","-"))
x
ranges(x)
plotGRanges = function(x) plotRanges(ranges(x))
par(mfrow=c(2,1))
plotGRanges(x)
plotGRanges(resize(x,1))

x = GRanges("chr1", IRanges(c(101,201,401,501),c(150,250,450,550)), strand="+")
y = GRanges("chr1", IRanges(c(101,221,301,401,541),c(150,250,350,470,550)), strand="+")
par(mfrow=c(2,1))
GRangesList(x,y)
plotGRanges(c(x,y))
c(x,y)
width(c(x,y))
xy <- disjoin(c(x,y))
sum(width(xy[xy %over% x]))
sum(width(xy[xy %over% y]))
210-140

disjoined = disjoin(c(x,y))

not.in.both = !(disjoined %over% x & disjoined %over% y)

sum(width(disjoined[ not.in.both ]))
z <- GRanges(x, strand = "-")
x %over% z


#2020-07-28 
library(ERBS)
data(HepG2)
data(GM12878)
HepG2[17]
index <- distanceToNearest(HepG2[17], GM12878)
GM12878[945]


res <- distanceToNearest(HepG2, GM12878)
d <- as.matrix(mcols(res))
l <- length(d)
l1 <- length(d[d<2000])
prop <-  l1/l
prop  
  
  
#2020-08-10
#I had to jump one section because the course expires for me 
library(BSgenome)
library(Biostrings)
ag = available.genomes()
ag 
grep("rerio", ag, value=TRUE)
grep("Hsap", ag, value=TRUE)

# inspect the human genome
library(BSgenome.Hsapiens.UCSC.hg19)
Hsapiens
length(Hsapiens)
class(Hsapiens)
methods(class="BSgenome")

# inspect human genome
Hsapiens$chrX
substr(Hsapiens$chrX, 5e6, 5.1e6)
nchar(Hsapiens$chrY)
nchar(Hsapiens[[24]])

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19.masked")
library(BSgenome.Hsapiens.UCSC.hg19.masked)
c17m <- BSgenome.Hsapiens.UCSC.hg19.masked$chr17
class(c17m)
c22m <- BSgenome.Hsapiens.UCSC.hg19.masked$chr22
c22m
3.198546e-01 * 100


# inspect data available from Ensembl
BiocManager::install("ensembldb")
library(ensembldb)
BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
names(listTables(EnsDb.Hsapiens.v75))

# extract Ensembl transcripts
edb = EnsDb.Hsapiens.v75  # abbreviate
txs <- transcripts(edb, filter = GeneNameFilter("ZBTB16"),
                   columns = c("protein_id", "uniprot_id", "tx_biotype"))
txs

# compare Ensembl and UCSC transcripts
alltx = transcripts(edb)    # Ensembl is larger
utx = transcripts(txdb)    # UCSC is smaller

# table of biological types of transcripts
table(alltx$tx_biotype)


BiocManager::install("Gviz")
library(ph525x)
stopifnot(packageVersion("ph525x") >= "0.0.16") # do over if fail 
modPlot("ESR1", useGeneSym=FALSE, collapse=FALSE) 

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

library("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
select(org.Hs.eg.db, keys="ESR1", keytype = "SYMBOL", columns ="ENTREZID")

length(transcripts(txdb,filter = list(gene_id = "2099" ))) 
