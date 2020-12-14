# Rpackages
source("http://bioconductor.org/biocLite.R")
if (!(require(annotate))) biocLite("annotate")
if (!(require(GOstats))) biocLite("GOstats")
if (!(require(org.Mm.eg.db))) biocLite("org.Mm.eg.db")
# Data
## Read data
topTab <- read.table("http://ueb.vhir.org/tiki-download_file.php?fileId=2689", 
                     head=TRUE, sep=",", dec=".")
expres <- read.table("http://ueb.vhir.org/tiki-download_file.php?fileId=2690", 
                     head=TRUE, sep="\t", dec=".")

## Unique identifiers
geneList <- as.character(topTab$SYMBOL); length(geneList)
backgrd <- as.character(expres$SYMBOL); length(backgrd)

## Symbols2Entrezs
require(org.Mm.eg.db)
geneListEntrezs <- unlist(mget(geneList, org.Mm.egSYMBOL2EG, ifnotfound=NA))
backgrdEntrezs <-  unlist(mget(backgrd, org.Mm.egSYMBOL2EG, ifnotfound=NA)) 

geneListEntrezs <- geneListEntrezs[!is.na(geneListEntrezs)]
backgrdEntrezs <- backgrdEntrezs[!is.na(backgrdEntrezs)]

# GOAnalysis
require(GOstats)
geneIds <- unique(geneListEntrezs)
entrezUniverse <- unique(backgrdEntrezs)

## Creamos los "hiperparametros" en que se basa el analisis
GOparams = new("GOHyperGParams",
               geneIds=geneIds, universeGeneIds=entrezUniverse,
               annotation="org.Mm.eg.db", ontology="BP",
               pvalueCutoff=0.001, conditional=FALSE,
               testDirection="over")

## Ejecutamos los analisis
  
GOhyper = hyperGTest(GOparams)
cat("GO\n")
print(head(summary(GOhyper)))

# Creamos un informe html con los resultados
GOfilename =file.path(paste("GOResults.",".html", sep=""))
htmlReport(GOhyper, file = GOfilename, summary.args=list("htmlLinks"=TRUE))

