# Visualize spatial data
# Following Seurat workflow: https://satijalab.org/seurat/archive/v3.2/spatial_vignette.html

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DropletQC)

# Reload
sobj <- readRDS("~/Dropbox/Zoe/scf_version/analysis/spatial/merged_spatial_sobj.RDS")

# Load data
# Note: need filtered matrices, h5, and spatial folder, I think
sobj.data <- Load10X_Spatial("spatial/L212_B1_VIS/outs/")

# Check QC
VlnPlot(sobj.data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
SpatialFeaturePlot(sobj.data, features = "nCount_Spatial") + theme(legend.position = "right")

# Use scTransform
sobj <- SCTransform(sobj.data, assay = "Spatial") %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:30) %>%
  RunTSNE()

sobj <- FindClusters(sobj, resolution = 0.4)

# Check out marker genes on spatial plot
sobj <- PercentageFeatureSet(sobj, pattern = "^MT-", col.name = "percent.mt")
sobj <- PercentageFeatureSet(sobj, pattern = "^RPL|^Rpl|^RPS|^Rps", col.name = "percent.ribo")
SpatialFeaturePlot(sobj, features = c("C8G"))
SpatialFeaturePlot(sobj, features = "nFeature_Spatial")
SpatialFeaturePlot(sobj, features = "percent.mt")

# View directly on raw slice
SpatialFeaturePlot(sobj, features = c("SERPINA1"),
                   pt.size.factor = 1, ncol = 1, crop = FALSE)
                   #alpha = c(0.1, 1))

# To calculate per-cluster expression stats
medEx <- as.data.frame(cbind(sobj$percent.ribo, Idents(sobj)))
colnames(medEx) <- c("expression","cluster")
medEx %>% group_by(cluster) %>% summarise(median(expression))

# Visualize as single cell plots
DimPlot(sobj, reduction = "umap", label = TRUE)
DimPlot(sobj, reduction = "tsne", label = TRUE)
FeaturePlot(sobj, features = "SERPINA1", reduction = "umap")
FeaturePlot(sobj, features = "nCount_SCT", reduction = "umap")

# Find markers - this will take a while
sobj <- FindSpatiallyVariableFeatures(sobj, assay = "SCT", features = VariableFeatures(sobj)[1:1000], 
                              selection.method = "markvariogram")

# Get top markers
top.features <- head(SpatiallyVariableFeatures(sobj, selection.method = "markvariogram"), 15)
top.features
# Find markers between clusters 5 and 6
de_markers <- FindMarkers(sobj, ident.1 = 0, ident.2 = 1, recorrect_umi = FALSE, only.pos = TRUE, 
                          min.pct = 0.1, logfc.threshold = 0.2)
# Find all markers
all.markers <- FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,
                              recorrect_umi = FALSE)
top.markers <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# Do other spatial stuff
SpatiallyVariableFeatures(sobj, selection.method = "markvariogram")[1:50]
SpatialFeaturePlot(sobj, features = top.features, ncol = 3, alpha = c(0.1, 1))

# Visualize clusters on spatial plot
SpatialDimPlot(sobj)

# Violin plots of markers per cluster
VlnPlot(sobj, features = "C8G")
DimPlot(sobj, reduction = "umap")
FeaturePlot(sobj, features = "C8G")

# Save
saveRDS(sobj, file = "~/Dropbox/Zoe/scf_version/analysis/spatial/L212_B1_VIS/B1_spatial_sobj.RDS")

# Reload
sobj1 <- readRDS(file = "~/Dropbox/Zoe/scf_version/analysis/spatial/L212_A1_VIS/A1_spatial_sobj.RDS")
sobj2 <- readRDS(file = "~/Dropbox/Zoe/scf_version/analysis/spatial/L212_B1_VIS/B1_spatial_sobj.RDS")
# Merge spatial objects together
# Follow vignette: https://satijalab.org/seurat/articles/spatial_vignette.html 
sobj <- merge(sobj1, sobj2)
DefaultAssay(sobj) <- "SCT"
VariableFeatures(sobj) <- c(VariableFeatures(sobj1), VariableFeatures(sobj2))
sobj <- RunPCA(sobj, verbose = FALSE)
sobj <- FindNeighbors(sobj, dims = 1:30)
sobj <- FindClusters(sobj, verbose = FALSE, resolution = 0.2)
sobj <- RunUMAP(sobj, dims = 1:30)
# Save
saveRDS(sobj, file = "~/Dropbox/Zoe/scf_version/analysis/spatial/merged_spatial_sobj.RDS")
sobj <- readRDS("~/Dropbox/Zoe/scf_version/analysis/spatial/merged_spatial_sobj.RDS")

dataset <- sobj
markers <- FindAllMarkers(dataset,
                          slot = "scale.data", assay = "SCT",
                          only.pos = TRUE, recorrect_umi = FALSE)
sigMarkers <- markers[markers$p_val_adj < 0.05,]
# Check how many markers are in dataset
for (clust in levels(as.factor(sigMarkers$cluster))) {
  print(paste("Cluster", clust, "has", sum(sigMarkers$cluster == clust), "DE genes."))
}
write.table(sigMarkers,
            file = "~/Dropbox/Zoe/scf_version/analysis/spatial/merged_spatial_res02_markers.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Get DE vs rest for all clusters
# Set identity
Idents(sobj) <- "SCT_snn_res.0.2"
# Get DE counts from clusters for pathway analysis
markers <- FindAllMarkers(sobj,
                        only.pos = FALSE,
                        min.pct = 0, 
                        logfc.threshold = 0,
                        assay = "SCT",
                        slot = "counts",
                        recorrect_umi = FALSE)
# OR do one clust vs another BUT DOESN'T REALLY WORK says me a long time ago
markers <- FindMarkers(sobj,
                       ident.1 = "0",
                       ident.2 = "1",
                       only.pos = FALSE,
                       min.pct = 0,
                       logfc.threshold = 0,
                       assay = "SCT",
                       slot = "counts",
                       recorrect_umi = FALSE)
# Get desired cluster (If used FindAllMarkers)
clust <- markers[markers$cluster == 0,]
# Keep only avg_diff and rename column
clust <- dplyr::select(clust, avg_log2FC)
# Add gene name as a column
clust$gene <- row.names(clust)
# Sort
clust <- clust[order(-clust$avg_log2FC),]
# Reorder columns
clust <- dplyr::select(clust, gene, avg_log2FC)
write.table(clust, file = "~/Dropbox/Zoe/scf_version/analysis/spatial/combined_pathways_res0.2/0_vs_1.rnk",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Now for correlation tests
# Just take averages across clusters below
Idents(sobj) <- "SCT_snn_res.0.2"
spatialAv <- AverageExpression(sobj, assays = "SCT", slot = "scale.data")$SCT
# Now compare to TLH samples
load("integration/L212_L202_TLH_rpca_mito50_scClustViz.RData")
Idents(scSeurat) <- "integrated_snn_res.0.6"
scAv <- AverageExpression(scSeurat, assays = "SCT", slot = "scale.data")$SCT
# Choose clusters (if desired)
spatialAv <- spatialAv[,c("0","1")]
# Scale and stuff
spatialAv <- as.data.frame(na.omit(t(scale(t(spatialAv)))))
scAv <- as.data.frame(na.omit(t(scale(t(scAv)))))
# Order by row name
spatialAvOrd <- spatialAv[order(row.names(spatialAv)),]
scAvOrd <- scAv[order(row.names(scAv)),]
# Now find intersecting genes
matches <- intersect(row.names(spatialAvOrd), row.names(scAvOrd))
# Look at how many genes matched
length(matches)
# Make new matrices with only matching gene names
spatialAvCor <- spatialAvOrd[matches,]
scAvCor <- scAvOrd[matches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(spatialAvCor))), seq(ncol(as.matrix(scAvCor)))),
            mapply(function(i, j) cor.test(as.matrix(spatialAvCor)[, i], as.matrix(scAvCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(spatialAvCor)),
               dimnames=list(colnames(as.matrix(spatialAvCor)), colnames(as.matrix(scAvCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(spatialAvCor)),
                dimnames=list(colnames(as.matrix(spatialAvCor)), colnames(as.matrix(scAvCor))))
heatmap(cors,
        main = "Pearson of A1 spatial slice vs integrated single-cell",
        xlab = "Single-cell integrated clusters",
        ylab = "Spatial A1 clusters")

# Deconvolution with integrated data
# Following https://satijalab.org/seurat/articles/spatial_vignette.html
sobj <- readRDS(file = "spatial/L212_A1_VIS/A1_spatial_sobj.RDS")
Idents(sobj) <- "SCT_snn_res.0.2"
load("L212/TLH/L212_TLH_cellranger_mito50_scClustViz.RData")
Idents(scSeurat) <- "SCT_snn_res.0.4"
# Take a look at the single-cell data
DimPlot(scSeurat, group.by = "SCT_snn_res.0.4", label = TRUE)
# Find anchors between single-cell and spatial data
anchors <- FindTransferAnchors(reference = scSeurat, query = sobj, normalization.method = "SCT")
# Transfer single-cell clusters to spatial
predictions.assay <- TransferData(anchorset = anchors, refdata = scSeurat$SCT_snn_res.0.4, prediction.assay = TRUE,
                                  weight.reduction = sobj[["pca"]], dims = 1:30)
# Add predictions and make default assay
sobj[["predictions"]] <- predictions.assay
DefaultAssay(sobj) <- "predictions"
# Visualize certain single-cell clusters as predicted on spatial data
SpatialFeaturePlot(sobj, features = "2")
# Or multiple at a time
SpatialFeaturePlot(sobj, features = c("0", "1", "3"), pt.size.factor = 1.6, ncol = 3, crop = TRUE)
# Now... do something different. Find top clusters.
sobj <- FindSpatiallyVariableFeatures(sobj, assay = "predictions", selection.method = "markvariogram",
                                        features = rownames(sobj), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(sobj), 4)
SpatialPlot(object = sobj, features = top.clusters, ncol = 2)
# Project onto liver slice
SpatialFeaturePlot(sobj, features = c("10","4","5","20"),
                   pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))

# Merge spatial objects together
# Follow vignette: https://satijalab.org/seurat/articles/spatial_vignette.html 
sobj1 <- readRDS("spatial/L212_A1_VIS/A1_spatial_sobj.RDS")
sobj.merge <- merge(sobj1, sobj2)

# Can I run dropletQC on spatial data?
nf1 <- nuclear_fraction_tags(
  outs = "~/Dropbox/Zoe/scf_version/analysis/spatial/L212_A1_VIS/outs",
  tiles = 20, cores = 1, verbose = FALSE)
# Add _1 to cell names
rownames(nf1) <- gsub("$","_1",rownames(nf1))
# Now do other slice
nf2 <- nuclear_fraction_tags(
  outs = "~/Dropbox/Zoe/scf_version/analysis/spatial/L212_B1_VIS/outs",
  tiles = 20, cores = 1, verbose = FALSE)
rownames(nf2) <- gsub("$","_2",rownames(nf2))
all_nf <- rbind(nf1,nf2)
# Subset nfs to get matching order
nuclear_fraction <- all_nf[row.names(sobj@meta.data),]
# Add this as metadata
sobj <- AddMetaData(sobj, nuclear_fraction, col.name = "nuclear_fraction")
# Visualize
SpatialFeaturePlot(sobj, features = "nuclear_fraction")
FeaturePlot(sobj, features = "nuclear_fraction")
VlnPlot(sobj, features = "nuclear_fraction")
FeatureScatter(sobj, feature1 = "nuclear_fraction", feature2 = "APOA1")
FeatureScatter(sobj, feature1 = "nuclear_fraction", feature2 = "CYP2E1")
FeatureScatter(sobj, feature1 = "nuclear_fraction", feature2 = "C8G")
FeatureScatter(sobj, feature1 = "nuclear_fraction", feature2 = "percent.mt")
# Pericentral heps have a slightly higher nuclear fraction; can kind of see
# with violin plot and by comparing nuclear fraction to periportal markers
# on spatial plot

