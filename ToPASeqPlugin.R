## ----setup, echo=FALSE--------------------------------------------------------
suppressPackageStartupMessages({ 
    library(ToPASeq)
    library(EnrichmentBrowser)
    library(graphite)
    library(BiocStyle)
})

## ----lib----------------------------------------------------------------------
library(ToPASeq)


input <- function(inputfile) {
## ----loadAirway---------------------------------------------------------------
#library(airway)
#data(airway)
mydata <<- readRDS(inputfile)
}

run <- function() {}

output <- function(outputfile) {
## ----processAirway------------------------------------------------------------
airSE <- mydata[grep("^ENSG", rownames(mydata)),]
dim(airSE)
assay(airSE)[1:4,1:4]

## ----pdataAirway--------------------------------------------------------------
airSE$GROUP <- ifelse(mydata$dex == "trt", 1, 0)
table(airSE$GROUP)

## ----pdataAirway2-------------------------------------------------------------
airSE$BLOCK <- mydata$cell
table(airSE$BLOCK)

## ----deAirway-----------------------------------------------------------------
library(EnrichmentBrowser)
airSE <- deAna(airSE, de.method="edgeR")
rowData(airSE, use.names=TRUE)

## ----pwys---------------------------------------------------------------------
library(graphite)
pwys <- pathways(species="hsapiens", database="kegg")
pwys

## ----nodes--------------------------------------------------------------------
nodes(pwys[[1]])

## ----mapIDs-------------------------------------------------------------------
airSE <- idMap(airSE, org="hsa", from="ENSEMBL", to="ENTREZID")

## ----genes--------------------------------------------------------------------
all <- names(airSE)
de.ind <- rowData(airSE)$ADJ.PVAL < 0.01
de <- rowData(airSE)$FC[de.ind]
names(de) <- all[de.ind]

## ----nrGenes------------------------------------------------------------------
length(all)
length(de)

## ----prs----------------------------------------------------------------------
res <- prs(de, all, pwys[1:100], nperm=100)
head(res)

## ----prsWeights---------------------------------------------------------------
ind <- grep("ErbB signaling pathway", names(pwys))
weights <- prsWeights(pwys[[ind]], de, all)
write.csv(weights, outputfile)

## ----maxWeight----------------------------------------------------------------
weights[weights == max(weights)]

}
