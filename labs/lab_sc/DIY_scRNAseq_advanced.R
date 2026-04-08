# Introduction ----
# the goal of this script is to explore several advanced topics in single cell analysis, including working with multi-modal data (in this case, CITE-seq), subclustering, pseudobulking, and GO analysis.
# To explore these topics, we need yet another example dataset.  
# We'll use unpublished data from my lab that draws on how the mucosal immune system responds to different kinds of infections.
# You can download this dataset from the course lecture page or the data page.

# Lecture 1 major topics Multimodal data, Subclustering, Pseudobulk DEG, GO analysis 
# Load the Libraries----
library(Seurat) # does all the heavy lifting in this script, and needed to work with Seurat Objects
library(tidyverse)
library(gprofiler2) # used to convert human cell cycle genes to mouse orthologs
library(scCustomize) # alternative and often improved visuals for Seurat graphing functions
library(presto) # Seurat implements presto to enhance speed of "FindMarkers" type functions
library(plotly) # we'll use it here to make interactive UMAPs
library(ggrepel) # helps with our plot labels on a crowded UMAP space


# Load the Seurat Object ----
# In this experiment we are looking at the mesenteric lymph node from mice that are controls or that have been infected by one of several pathogens/microbes. For expediency, the cells in this project have already have preliminary coarse labels and fine/granular cell type labels. Additionally, the total cells have been dramatically reduced to work within the confines of a single class and resource limitations. 
# MLN_seurat <- LoadSeuratRds("MLN_subsampled_DIY.Rds")
lleum_seurat <- LoadSeuratRds("~/MBMI/CAMB7140/labs/lab_sc/Ileum_lab14.Rds")

# Now we normalize the Seurat object we just created and find variable features
lleum_seurat <-NormalizeData(lleum_seurat, verbose = FALSE)
lleum_seurat <- FindVariableFeatures(lleum_seurat, verbose = FALSE)

# Let's calculate percent of mitochondrial reads
# NOTE: change 'MT' to 'mt' for mouse
lleum_seurat[["percent.mt"]] <- PercentageFeatureSet(object = lleum_seurat, pattern = "mt") 
# in the violin plot above, features = genes detected, while counts = total molecules detected
# Make violin plot
VlnPlot(lleum_seurat, c("nCount_RNA", "nFeature_RNA", "percent.mt"), pt.size = 0.1)
# Filter your data
lleum_seurat <- subset(lleum_seurat, subset = 
                           nCount_RNA < 20000 & 
                           nCount_RNA > 1000 & 
                           nFeature_RNA > 1000 & 
                           percent.mt < 20)

# another QA plot
FeatureScatter(lleum_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# Potential things to look for in the type of QA plot produced above:
# 1. Data points in the bottom LEFT hand quadrant = low genes and UMIs per cell. May represent poor quality cells.
# 2. Data points in the bottom RIGHT hand quadrant = low genes but high UMIs per cell. These could be dying cells, but also could represent a population of a low complexity celltype (i.e red blood cells).

# take a look at the top variable features 
top10 <- head(VariableFeatures(lleum_seurat), 10)
top10

# Plot UMAP ----
# it is standard practice to apply a linear transformation ('scaling') before PCA. For single cell data this includes:
# 1. Shifting the expression of each gene, so that the mean expression across cells is 0
# 2. Scaling the expression of each gene, so that the variance across cells is 1
# This gives equal weight in downstream analyses, so that highly-expressed genes do not dominate
all.genes <- rownames(lleum_seurat)
lleum_seurat <- ScaleData(lleum_seurat, features = all.genes)
lleum_seurat <- RunPCA(lleum_seurat, npcs = 40, verbose = FALSE)
# What contributes to the PCA?  Let's take a look using the VizDimLoadings function
VizDimLoadings(lleum_seurat, dims = 1:2, reduction = "pca")
# How many dimensions should we keep? Usually not all 40 
ElbowPlot(lleum_seurat) #in this case I wouldn't keep more than 20 for the steps below

lleum_seurat <- RunUMAP(lleum_seurat, reduction = "pca", dims = 1:20)
lleum_seurat <- FindNeighbors(lleum_seurat, reduction = "pca", dims = 1:20)
lleum_seurat <- FindClusters(lleum_seurat, resolution = 0.5)
DimPlot(lleum_seurat, reduction = "umap", split.by = "orig.ident", label = TRUE)


# Find cluster-specific genes ----
# generally speaking there are three main ways you can find cluster-specific marker genes with Seurat
# 1. 'FindMarkers' to compare a select cluster to all other cells not in that cluster
# 2. 'FindAllMarkers' to compare EACH cluster to all other cells
# 3. 'FindConservedMarkers' to identify genes conserved (shared) between two defined clusters

# We'll start with FindMarkers, since it allows you to choose exactly which cluster you'd like to focus on.
# Seurat has implemented a speed-saving measure using presto package. Have them install this package and load in the beginning
library(presto) # install with devtools::install_github("immunogenomics/presto")
cluster1.markers <- FindMarkers(lleum_seurat, ident.1 = 1, min.pct = 0.25)
cluster1.markers$pct.diff <- cluster1.markers$pct.1 - cluster1.markers$pct.2
cluster1.markers.df <- as_tibble(cluster1.markers, rownames = "geneID")
# Export DEGs for each cluster (ranked by avg_logFC > 0.5)
myTopHits_cluster1 <- cluster1.markers.df %>% arrange(desc(avg_log2FC))
myTopHits_cluster1 <- dplyr::slice(myTopHits_cluster1, 1:20)

# you can make this list of genes into an interactive table
datatable(myTopHits_cluster1, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Cluster 1 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

# plot genes of interest on UMAP
FeaturePlot(lleum_seurat, 
            reduction = "umap", 
            features = c("Lef1"),
            pt.size = 0.4, 
            order = TRUE,
            #split.by = "orig.ident",
            min.cutoff = 'q10',
            label = FALSE)

# now let's try with FindAllMarkers
lleum.markers <- FindAllMarkers(lleum_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# let's take the top 10 marker genes for each cluster and plot as a heatmap
top10 <- lleum.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

DoHeatmap(lleum_seurat, features = top10$gene)

# Assigning identity to cell clusters  ----
library(SingleR) #automated cell type annotation ('label transfer') using reference data
library(celldex) #a large collection of reference expression datasets with curated cell type labels for use with SingleR package
library(pheatmap)
library(Azimuth) # providing the code below for help with Azimuth installation - Say "yes" to source installation
# remotes::install_github('satijalab/azimuth', ref = 'master')
# depending on your existing package environment, unfortunately this can sometimes be a headache to install
# If azimuth fails to install, you can try the following
# set a higher timeout option for R as follows:
# options(timeout = max(1000, getOption("timeout")))
# this high timeout will let you install the following dependency, if R is complaining that it wants it
# install.packages("BSgenome.Hsapiens.UCSC.hg38")
# now try to install Azimuth again
# remotes::install_github('satijalab/azimuth', ref = 'master')
# if running Azimuth fails and complains with the error "object 'CRsparse_colSums' not found", this is a documented error (https://github.com/satijalab/seurat/issues/8202#issue-2047511055)
# you can fix this with by reinstalling TFBStools as follows:
# BiocManager::install("TFBSTools", type = "source", force = TRUE)

# it can also be useful to turn the Seurat object into a singleCellExperiment object, for better interoperability with other bioconductor tools
# two ways to get singleCellExperiment object

# option 2 - Seurat allows you to convert directly
lleum.sce <- as.SingleCellExperiment(lleum_seurat)

# the singleCellExperiment data structure is easy to work with
rownames(lleum.sce)
colnames(lleum.sce)
reducedDims(lleum.sce)
assays(lleum.sce)
my.subset <- lleum.sce[,c(1,2,8)]
rowData(lleum.sce)$Symbol <- rownames(lleum.sce)

# OPTION 1: assign identity to cell clusters using public RNA-seq datasets
# To do this, we'll use singleR and celldex (requires an internet connection to connect to ExperimentHub)
ENCODE.data <- BlueprintEncodeData(ensembl = FALSE) #259 RNA-seq samples of pure stroma and immune cells as generated and supplied by Blueprint and ENCODE
HPCA.data <- HumanPrimaryCellAtlasData(ensembl = FALSE) #713 microarray samples from the Human Primary Cell Atlas (HPCA) (Mabbott et al., 2013).
DICE.data <- DatabaseImmuneCellExpressionData(ensembl = FALSE) #1561 bulk RNA-seq samples of sorted immune cell populations
ImmGen.data <- ImmGenData(ensembl = FALSE) # 830 microarray samples of pure mouse immune cells, generated by the Immunologic Genome Project (ImmGen)
Monaco.data <- MonacoImmuneData(ensembl = FALSE) #114 bulk RNA-seq samples of sorted immune cell populations that can be found in GSE107011.
MouseRNAseq.data <- MouseRNAseqData(ensembl = FALSE) #358 bulk RNA-seq samples of sorted cell populations
Hemato.data <- NovershternHematopoieticData(ensembl = FALSE) #211 bulk human microarray samples of sorted hematopoietic cell populations that can be found in GSE24759

predictions <- SingleR(test=lleum.sce, assay.type.test=1, 
                       ref=Monaco.data, labels=Monaco.data$label.main)

plotScoreHeatmap(predictions)

#now add back to singleCellExperiment object (or Seurat objects)
lleum.sce[["SingleR.labels"]] <- predictions$labels
plotUMAP(lleum.sce, colour_by = "SingleR.labels")