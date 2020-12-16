
##Load the data
library("airway")

#to find out where on your computer the files from a package have been installed.
indir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(indir)

#read the targets file of the experiment
csvfile <- file.path(indir, "sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
sampleTable

##store the bam files in object filenames
filenames <- file.path(indir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)
filenames

##Indicate that these are bam files
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)


##Defining gene models
library("GenomicFeatures")

gtffile <- file.path(indir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format = "gtf", circ_seqs = character())
txdb

##resume the exons by gene
ebg <- exonsBy(txdb, by="gene")
ebg

##Generating count Matrices
library("GenomicAlignments")

se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )

##See the count matrix
se
assay(se)

#dimensions of the count matrix
dim(se)

#############################################################################
##START FROM THE COMPLETE DATA
#############################################################################
BiocManager::install("airway")
library(airway)
data("airway")
se <- airway
dim(se)


#reorder the levels 
library(magrittr)
se$dex %<>% relevel("untrt")
se$dex

##Milions of fragments aligned to the genes
round( colSums(assay(se)) / 1e6, 1 )

##Check that the object is correct
colData(se)

##Make a DESeqDataSet object
library("DESeq2")

dds <- DESeqDataSet(se, design = ~ cell + dex)


##Filter out those rows without any count
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

##variance stabilizing transformation?(VST)
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)


##Sample distances.
sampleDists <- dist(t(assay(vsd)))
sampleDists

##Heatmap + Hierarchical clustering
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$dex, vsd$cell, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

##Principal Component Analysis Plot
plotPCA(vsd, intgroup = c("dex", "cell"))


##Differential Expression Analysis
dds <- DESeq(dds)

##building results table
res <- results(dds, contrast=c("dex","trt","untrt"))
##Toptable
head(res[order(res$padj),])

##Information about the columns of the results table
mcols(res, use.names = TRUE)

##Subset the significant genes with strongest downregulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

##Subset the significant genes with strongest upregulation
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])


##Heatmap of genes
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)

mat  <- assay(vsd)[topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
pheatmap(mat, annotation_col = anno)


##Results annotation
library("org.Hs.eg.db")
library("AnnotationDbi")

res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")


resOrdered <- res[order(res$pvalue),]
head(resOrdered)

##Exporting the results
resOrderedDF <- as.data.frame(resOrdered)[1:100, ]
write.csv(resOrderedDF, file = "results.csv")
