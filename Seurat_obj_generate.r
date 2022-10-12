############################### Aggregating data after STARSolo and generating Seurat object
#########################
#########################
######### snRNA Multiome samples 

### request small amount of resources
srun --mem=32GB --time=6:00:00 --pty --cpus-per-task=4 --x11 -p bigmem,interactive bash
## activate conda environment
module load miniconda/4.12.0
# source activate /gpfs/gibbs/pi/kaminski/ab3478/conda_envs/agshin_R4env
source activate agshin_R4env

R

## load conda to enter R environment
library(Matrix)
library(tidyverse)
library(ggplot2)
library(scales)
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(viridisLite)
library(circlize)
options(stringsAsFactor=FALSE)


########################## knee plot function
betterBarcode_rank_plot <- function(sce, nBarcodes, sample){
    counts <- sce$metadata$nUMI[order(sce$metadata$nUMI, decreasing=TRUE)]
    counts <- sort(counts, decreasing = TRUE)
    ranks <- seq(length(counts))
    datf <- data.frame(Rank = ranks, Count = counts)
    datf <- datf[!is.na(datf[, 2]), , drop = FALSE]

    counts[duplicated(counts)] <- NA
    ranks[duplicated(counts)] <- NA
    datf.plot <- data.frame(Rank = ranks, Count = counts)
    datf.plot <- datf.plot[!is.na(datf.plot[, 2]), , drop = FALSE]
    png(paste("Elbow.nUMI.total.", sce$metadata$library.ident, ".png", sep=""), res=200, unit="in", height=8, width=11)
    print({
        p <- ggplot(datf.plot, aes_string(x = "Rank", y = "Count")) +
        geom_point() +
        scale_x_continuous(trans = "log10", labels = comma,
            breaks=c(10, 100, 500, 1000, 5000, 10000, 30000, 50000, 100000, 500000),
            limits=c(1, 500000)) +
        scale_y_continuous(name = "Droplet size", trans = "log10",labels = comma,
            breaks=c(1, 10, 100, 150, 200, 500, 1000, 5000, 10000, 25000, 50000),
            limits=c(1, 50000)) +
            geom_vline(xintercept=nBarcodes, linetype="dashed", color="red")  +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(sample)
    })
    dev.off()
    return(datf)
}


## this the parent directory for the sample outputs
topLevel.dir <- "/vast/palmer/scratch/kaminski/ab3478/snRNA_Multiome/sample_out"
### we will work in analysis folder
dir.create("/vast/palmer/scratch/kaminski/ab3478/snRNA_Multiome/analysis")
setwd("/vast/palmer/scratch/kaminski/ab3478/snRNA_Multiome/analysis")

sample.paths <- list.files(topLevel.dir, full.names=TRUE)
sample.names <- basename(sample.paths)

##### loop over samples, pull the filtered matrix, append sample identifier, collapse the matrices
### make an empty list
list.counts <- list()
list.metadata <- list()

collectSTARsolo <- function(sample.path, nBarcodes){
    sample.name <- basename(sample.path)
    cat("collecting \'GeneFull\' count information\n")
#### get the gene-full information first
    geneFull.path <- file.path(sample.path, "Solo.out", "GeneFull","raw")
    gf.temp.barcodes <- read.table(file.path(geneFull.path, "barcodes.tsv"), sep="\t", header=FALSE)[,1]
    gf.temp.genes <- as.vector(read.table(file.path(geneFull.path, "features.tsv"), sep="\t", header=FALSE)$V2)
    gf.counts <-  Matrix::readMM(file.path(geneFull.path, "matrix.mtx"))
### put them together
    rownames(gf.counts) <- gf.temp.genes
    colnames(gf.counts) <- gf.temp.barcodes
############################################
########### gathering the rest of the stuff
    cat("collecting \'Gene\' and \'Velocyto\' splice information\n")
    gene.path <- file.path(sample.path, "Solo.out", "Gene","raw")
    g.temp.barcodes <- read.table(file.path(gene.path, "barcodes.tsv"), sep="\t", header=FALSE)[,1]
    g.temp.genes <- as.vector(read.table(file.path(gene.path, "features.tsv"), sep="\t", header=FALSE)$V2)
    g.counts <-  Matrix::readMM(file.path(gene.path, "matrix.mtx"))
### put them together
    rownames(g.counts) <- g.temp.genes
    colnames(g.counts) <- g.temp.barcodes
#### get % umi in Gene vs GeneFull
    percSpliced.g <- 100*(Matrix::colSums(g.counts)/Matrix::colSums(gf.counts))
## replace NA w/ zero
    percSpliced.g[is.nan(percSpliced.g)] <- 0
## replace Inf w/ 100
    percSpliced.g[is.infinite(percSpliced.g)] <- 100
################## get the scVelo spliced/unspliced estimate matrix
    velocyto.path <- file.path(sample.path, "Solo.out", "Velocyto","raw")
    velo.t <- readr::read_delim(file.path(velocyto.path, "matrix.mtx"), col_names = FALSE, delim = ' ', skip = 3)
    v.temp.barcodes <- read.table(file.path(velocyto.path, "barcodes.tsv"), sep="\t", header=FALSE)[,1]
#### split the three levels of information into separate matrices
    spliced.mtx <- Matrix::sparseMatrix(i = velo.t$X1, j = velo.t$X2, x = velo.t$X3, dims = dim(gf.counts))
    unspliced.mtx <- Matrix::sparseMatrix(i = velo.t$X1, j = velo.t$X2, x = velo.t$X4, dims = dim(gf.counts))
    ambiguous.mtx <- Matrix::sparseMatrix(i = velo.t$X1, j = velo.t$X2, x = velo.t$X5, dims = dim(gf.counts))
#### get an aggregated colsum of all barcodes
    all.colsums <- Matrix::colSums(spliced.mtx) + Matrix::colSums(unspliced.mtx) + Matrix::colSums(ambiguous.mtx)
#### assign the barcode names to the values
    percSpliced.v <- 100*(Matrix::colSums(spliced.mtx)/all.colsums)
    names(percSpliced.v) <- v.temp.barcodes
## replace NA w/ zero and Inf w/ 100, assign names
    percSpliced.v[is.nan(percSpliced.v)] <- 0
    percSpliced.v[is.infinite(percSpliced.v)] <- 100
    percUnSpliced.v <- 100*(Matrix::colSums(unspliced.mtx)/all.colsums)
    names(percUnSpliced.v) <- v.temp.barcodes
## replace NA w/ zero and Inf w/ 100, assign names
    percUnSpliced.v[is.nan(percUnSpliced.v)] <- 0
    percUnSpliced.v[is.infinite(percUnSpliced.v)] <- 100
## replace NA w/ zero and Inf w/ 100, assign names
    percAmbiguous.v <- 100*(Matrix::colSums(ambiguous.mtx)/all.colsums)
    percAmbiguous.v[is.nan(percAmbiguous.v)] <- 0
    percAmbiguous.v[is.infinite(percAmbiguous.v)] <- 100
    names(percAmbiguous.v) <- v.temp.barcodes
#### now only keep the barcodes selected from nBarcodes
#     percSpliced.v <- percSpliced.v[barcodes.keep]
#     percUnSpliced.v <- percUnSpliced.v[barcodes.keep]
#     percAmbiguous.v <- percAmbiguous.v[barcodes.keep]
#### collect the remaining metadata values of interest
    gf.nUMI <- Matrix::colSums(gf.counts)
    gf.nGenes <- Matrix::colSums(gf.counts > 0)
    mito.genes <- grep("^MT-", rownames(gf.counts), value=TRUE)
    percMito <- 100*(Matrix::colSums(gf.counts[mito.genes,])/gf.nUMI)
    percMito[is.nan(percMito)] <- 0
    percMito[is.infinite(percMito)] <- 100
    percMALAT1 <- 100*(gf.counts["MALAT1",]/gf.nUMI)
    percMALAT1[is.nan(percMALAT1)] <- 0
    percMALAT1[is.infinite(percMALAT1)] <- 100
### output each iteration into spot on list
    metadata <- data.frame(barcode=colnames(gf.counts),
        nUMI=gf.nUMI[colnames(gf.counts)], nGene=gf.nGenes[colnames(gf.counts)],
        percMito=percMito[colnames(gf.counts)], percMALAT1=percMALAT1[colnames(gf.counts)],
        percSpliced.g=percSpliced.g, percSpliced.v=percSpliced.v[colnames(gf.counts)],
        percUnspliced.v=percUnSpliced.v[colnames(gf.counts)], percAmbiguous.v=percAmbiguous.v[colnames(gf.counts)],
        library.ident=rep.int(sample.name, times=length(colnames(gf.counts))))
    out.list <- list(gf.counts, metadata)
    names(out.list) <- c("counts", "metadata")
    rank.matrix <- betterBarcode_rank_plot(out.list, nBarcodes, out.list$metadata$library.ident)
    umi.threshold <- rank.matrix[rank.matrix$Rank == nBarcodes,]$Count
    cat("nUMI threshold for library ", sample.names[i] , " = ", umi.threshold, "\n")
    rownames(out.list$metadata) <- paste(out.list$metadata$library.ident, "__", rownames(out.list$metadata),sep="")
    colnames(out.list$counts) <- rownames(out.list$metadata)
    out.list$metadata=out.list$metadata[out.list$metadata$nUMI >= umi.threshold,]
    cat("Number of barcodes extracted from library ", sample.names[i] , " = ", nrow(out.list$metadata), "\n")
    out.list$counts=out.list$counts[,colnames(out.list$counts) %in% rownames(out.list$metadata)]
    return(out.list)
}

for (i in 1:length(sample.paths)) {
  output <- collectSTARsolo(sample.paths[i], nBarcodes = 30000)
  list.counts[[i]] <- output$counts
  list.metadata[[i]] <- output$metadata
}

##Create big seurat Object
count.matrix=do.call(cbind, list.counts)
metadata.matrix = do.call(rbind, list.metadata)

dim(count.matrix)
# [1]  60651 121465
dim(metadata.matrix)
# [1] 121465     10
seurat_obj <- CreateSeuratObject(counts=count.matrix, meta.data=metadata.matrix, names.delim="__")

##Save seurat object
saveRDS(seurat_obj, file = "Seurat.obj.aggregate.snRNA.Multiome.101222.Rds")
