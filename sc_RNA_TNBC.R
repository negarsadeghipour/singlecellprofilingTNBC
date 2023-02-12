# In this script the single cell RNA-Seq results from this paper is regenerated:
# https://pubmed.ncbi.nlm.nih.gov/34956864/
# Bulk and Single-Cell Profiling of Breast Tumors Identifies TREM-1 as a 
#Dominant Immune Suppressive Marker Associated With Poor Outcomes


# data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188600
# using Seurat tutorial

#setwd("~/Documents/Bioinformatics pipelines with examples and R codes/Single Cell RNASeq/Projects/TNBC-human/single_cell_TNBC")


# 1. library loading ------------------------------------------------------

library(Seurat)
library(tidyverse)
library(GEOquery)



# 2. Fetch data from GEO --------------------------------------------------

# No need to perform the following two lines next time, as the folders are already created.
# getGEOSuppFiles("GSE188600")
# untar("GSE188600_RAW.tar", exdir = 'data/')

gds <- getGEO("GSE188600", GSEMatrix=TRUE)

# find the files in data folder (i.e., barcodes, features, matrix)
wd = getwd()
files = list.files(path = paste0(wd, '/data'), full.names = FALSE, recursive = FALSE)

# read the counts into a matrix
mtx.cnts <- ReadMtx(mtx = paste0('data/', files[3]),
                features = paste0('data/', files[2]),
                cells = paste0('data/', files[1]))
# create a seurat object
Seurat.obj <- CreateSeuratObject(counts = mtx.cnts)

Seurat.obj


# 3. QC and filtering -----------------------------------------------------

view(Seurat.obj@meta.data)

# calculate mitochondrial percentage (high percentage shows bad quality)
Seurat.obj$mito.prct <- PercentageFeatureSet(Seurat.obj, pattern = '^MT-')

# explore QC by two plots

VlnPlot(Seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "mito.prct"), ncol = 3)
FeatureScatter(Seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# filtering the cells to those with more than 800 RNA count and 500 genes and
# mitochondrial percentage less than 10
Seurat.obj.filt <- subset(Seurat.obj, subset = nCount_RNA >800 &
                          nFeature_RNA > 500 &
                            mito.prct < 10)



# 4. Finding Variable genes -----------------------------------------------

# Perform standard workflow steps to figure out if we see any batch effects
Seurat.obj.filt <- NormalizeData(object = Seurat.obj.filt)

# variable genes- with visualization of top 10 variable genes
Seurat.obj.filt <- FindVariableFeatures(object = Seurat.obj.filt, selection.method = 'vst', nfeatures = 2000)

top10 <- head(VariableFeatures(Seurat.obj.filt), 10)

plot1 <- VariableFeaturePlot(Seurat.obj.filt)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5. Clustering the cells -------------------------------------------------

# scaling the data
Seurat.obj.filt <- ScaleData(object = Seurat.obj.filt)
# perform linear dimensionality reduction
Seurat.obj.filt <- RunPCA(object = Seurat.obj.filt, features = VariableFeatures((Seurat.obj.filt)))
# selecting the PCA plots with elbow plot
ElbowPlot(Seurat.obj.filt, ndims = 35)
# clustering
Seurat.obj.filt <- FindNeighbors(object = Seurat.obj.filt, dim = 1:20)

# understanding the resolution
Seurat.obj.filt <- RunUMAP(object = Seurat.obj.filt, dim = 1:20)
Seurat.obj.filt <- FindClusters(object = Seurat.obj.filt, resolution = c(0.01, 0.1, 0.3, 0.5, 0.8, 1, 1.2))
DimPlot(Seurat.obj.filt, group.by = 'RNA_snn_res.0.5', label = TRUE)

Seurat.obj.filt <- RunUMAP(object = Seurat.obj.filt, dim = 1:20)
Seurat.obj.filt <- RunTSNE(object = Seurat.obj.filt, dims = 1:20)

# setting identity of clusters
Idents(Seurat.obj.filt)
Idents(Seurat.obj.filt) <- 'RNA_snn_res.0.5'

# plot to test for batch effects and finding optimal number of clusters
# If batch effects observed, integration is needed --- not in this example
# umap observatio
DimPlot(Seurat.obj.filt, reduction = 'umap', label = TRUE)
# tsne observation
DimPlot(Seurat.obj.filt, reduction = 'tsne' , label = TRUE)

# Find optimal number of clusters
DefaultAssay(Seurat.obj.filt) # make sure it is RNA

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Seurat.obj.filt.markers <- FindAllMarkers(Seurat.obj.filt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# select top (4 -> can be changed) upregulated genes in each cluster
clust.markers <- Seurat.obj.filt.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)



# 6. Steps to complete the annotation -------------------------------------

# top upregulated genes in each cluster
DoHeatmap(Seurat.obj.filt, features = clust.markers$gene, size = 4,
          angle = 90)
# DoHeatmap now shows a grouping bar, splitting the heatmap into groups or clusters. This can
# be changed with the `group.by` parameter
# total DE upregulated genes
DoHeatmap(Seurat.obj.filt, features = VariableFeatures(Seurat.obj.filt)[1:50], cells = 1:700, size = 2,
          angle = 90)


# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(Seurat.obj.filt, features = clust.markers$gene, ncol = 6)
# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(Seurat.obj.filt, features = clust.markers$gene, ncol = 4)
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(Seurat.obj.filt, features = clust.markers$gene, ncol = 4)
# Check individual feature
FeaturePlot(Seurat.obj.filt, features = c('CD163'), min.cutoff = 'q10')

# feature plots (feature identification based on up-reg genes and pangloadb)
features <- c("CST3", "CD86", "SPP1", "C1QA", "CTHRC1", "MGP", "CXCL13", "TNFRSF18", "JCHAIN", "IGKC",
              "NKG7", "GNLY", "IL7R", "RPL34", "STMN1", "KIAA0101")



# assigning the new clusters to data
new.cluster.ids <- c("Mature DC", "DC", "Fibroblast", "T cell", "Immature DC",
                     "NK cells", "Naive T cell", "TNBC")


names(new.cluster.ids) <- levels(Seurat.obj.filt)
Seurat.obj.filt <- RenameIdents(Seurat.obj.filt, new.cluster.ids)
DimPlot(Seurat.obj.filt, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(Seurat.obj.filt, reduction = "umap")

# Single cell heatmap of feature expression
DoHeatmap(subset(Seurat.obj.filt, downsample = 700), features = features, size = 3)
# Feature plot - visualize feature expression in low-dimensional space
FeaturePlot(Seurat.obj.filt, features = features, ncol = 4)

# select top (8 -> can be changed) upregulated genes in each cluster
top.features <- Seurat.obj.filt.markers %>%
  group_by(cluster) %>%
  slice_max(n = 6, order_by = avg_log2FC)
DoHeatmap(subset(Seurat.obj.filt, downsample = 700), features = top.features$gene, size = 3)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level
DotPlot(Seurat.obj.filt, features = features) + RotatedAxis()

# Gene Ontology -----------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library("EnsDb.Hsapiens.v86")

genes_to_test <- top.features$gene
mapIds <- mapIds(EnsDb.Hsapiens.v86, keys = genes_to_test, keytype = "SYMBOL", columns = c("GENEID"))
egIDs <- unlist(mget(mapIds, org.Hs.egENSEMBL2EG, ifnotfound = NA))

# KEGG
# converting to gene ID
kegg_results <- enrichKEGG(gene = egIDs, organism = "hsa", keyType = "kegg")
kegg_res_df <- as.data.frame(kegg_results)
fit_kegg <- plot(barplot(kegg_results, showCategory = 20))
dotplot(kegg_results)

# Gene Ontology 

GO_results <- enrichGO(gene = mapIds, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")
GO_res_df <- as.data.frame(GO_results)
dotplot(GO_results)
# x axis shows how many genes were involved in that pathway from the input data
fit_GO <- plot(barplot(GO_results, showCategory = 20))
dotplot(GO_results, showCategory = 20)

# Reactome

library(ReactomePA)
Reac_results <- enrichPathway(gene = egIDs, organism = "human", pAdjustMethod = "BH")
fit_Reac <- plot(barplot(Reac_results))
dotplot(Reac_results, showCategory = 20)

library(DOSE)
library(enrichplot)
edo <- enrichDGN(egIDs)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

library(msigdbr)
all_gene_sets <- msigdbr(species = "Homo sapiens")
msigdbr_df <- msigdbr(species = "human", category = "H")
msigdbr_t2g = msigdbr_df %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
test <- enricher(gene = genes_to_test, TERM2GENE = msigdbr_t2g)

# volcano plots
library(EnhancedVolcano)
Seurat.obj.filt.markers.posneg <- FindAllMarkers(Seurat.obj.filt, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

EnhancedVolcano(Seurat.obj.filt.markers.posneg,
                lab = rownames(Seurat.obj.filt.markers.posneg),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'Total differentially expressed genes',
                pCutoff = 10e-16,
                FCcutoff = 1.5,
                pointSize = 1,
                labSize = 3,
                col=c('gray', 'gray', 'gray', 'red3', 'blue'),
                colAlpha = 1,
                legendPosition = 'right')

# saving data
saveRDS(pbmc, file = "../output/pbmc3k_final.rds")

sessionInfo()
