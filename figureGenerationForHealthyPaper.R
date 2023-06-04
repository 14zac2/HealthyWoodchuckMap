# Load libraries

library(Seurat)
library(scClustViz)
library(ggplot2)
library(dplyr)

# Load desired datasets

# Liver map
sobj <- readRDS("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/allIntegrated_kanchor5_noBiopsyHeps_dropletQCFiltered.RDS")
Idents(sobj) <- "integrated_snn_res.0.4" # Lower resolution
Idents(sobj) <- "integrated_snn_res.1.2" # Higher resolution

# CCA Liver map
sobj <- readRDS("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/allIntegrated_cca_kanchor5_noBiopsyHeps_dropletQCFiltered.RDS")
Idents(sobj) <- "integrated_snn_res.1.4"

# Subclustered immune from liver
sobj <- readRDS("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/subclustering/immune_merged.RDS")
Idents(sobj) <- "SCT_snn_res.0.6"

# PBMCs
load("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/integrated_PBMC_kanchor5_scClustViz.RData")
sobj <- scSeurat
Idents(sobj) <- "integrated_snn_res.0.2"
Idents(sobj) <- "integrated_snn_res.0.8"

# CCA PBMCs
load("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/integrated_PBMC_cca_kanchor5_scClustViz.RData")
sobj <- scSeurat
Idents(sobj) <- "integrated_snn_res.0.6"

# Get some dotplots with specific markers

DotPlot(sobj,
        assay = "SCT",
        features = c("PTPRC", "CALCRL", "NKG7", "CD3E", "MARCO", "LYZ-1", "CD19", "MS4A1", "STAB2")
        ) +
  ggtitle("Select features for liver map")

# Get some feature plots for specific markers
mapType <- "Liver"

FeaturePlot(sobj, features = "sct_PTPRC") +
  ggtitle(paste("PTPRC -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CALCRL") +
  ggtitle(paste("CALCRL -", mapType, "map"))
FeaturePlot(sobj, features = "sct_NKG7") +
  ggtitle(paste("NKG7 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CD3E") +
  ggtitle(paste("CD3E -", mapType, "map"))
FeaturePlot(sobj, features = "sct_MARCO") +
  ggtitle(paste("MARCO -", mapType, "map"))
FeaturePlot(sobj, features = "sct_LYZ-1") +
  ggtitle(paste("LYZ -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CD19") +
  ggtitle(paste("CD19 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_MS4A1") +
  ggtitle(paste("MS4A1 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_STAB2") +
  ggtitle(paste("STAB2 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_ALB-1") +
  ggtitle(paste("ALB1 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CD4") +
  ggtitle(paste("CD4 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CD8A") +
  ggtitle(paste("CD8A -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CLEC4G") +
  ggtitle(paste("CLEC4G -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CD5L") +
  ggtitle(paste("CD5L -", mapType, "map"))
FeaturePlot(sobj, features = "sct_C1QB") +
  ggtitle(paste("C1QB -", mapType, "map"))
FeaturePlot(sobj, features = "sct_ACTA2") +
  ggtitle(paste("ACTA2 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_VWF") +
  ggtitle(paste("VWF -", mapType, "map"))
FeaturePlot(sobj, features = "sct_IGLL5-1") +
  ggtitle(paste("IGLL5 -", mapType, "map"))
FeaturePlot(sobj, features = "sct_CD68") +
  ggtitle(paste("CD68 -", mapType, "map"))
FeaturePlot(sobj, features = "LEF1") +
  ggtitle(expression(italic("LEF1")))
FeaturePlot(sobj, features = "RHEX") +
  ggtitle(expression(italic("RHEX")))

# Map with no cluster numbers
DimPlot(sobj, label = FALSE) & NoLegend()

# Map with cluster numbers
DimPlot(sobj, label = TRUE)

# Map with cell-type labels
DimPlot(sobj, group.by = "scina_labels_refined", label = TRUE) & NoLegend()

# Map by orig.ident
DimPlot(sobj, group.by = "orig.ident")

# Heatmap of main markers

# First get desired genes for liver map
plotFeatures <- c("CD3E", "FYN", "GZMK", "NKG7", # Lymphocyte
                  "IL7R-1", "IFITM3;IFITM2;IFITM1-2", # Lymphocyte
                  "BCL2","CD74", "VCAN", "XCL1;XCL2", "VWF", # Lymphocyte
                  "STAB2", "DNASE1L3", "CTSV;CTSL", "ITGA1", # Endothelial
                  "CD5L", "C1QC", "MSR1", "Clec4f", "S100A8", "S100A9", # Myeloid
                  "LYZ-1", "C1QB", "IL1B", # Myeloid
                  "MZB1", "JCHAIN", "IGLL5-1", # Antibody-secreting B cells
                  "ACTA2", "IGFBP7", "DCN", # Mesenchymal
                  "MS4A1", "CD79B", "BANK1", # B cells
                  "CFTR", "SLC4A4", "KRT19", # Cholangiocytes
                  "MKI67", "TOP2A", "STMN1", # Proliferative
                  "ALB-1", "APOC1", "APOE", # Hepatocytes
                  "CYP2E1", "FETUB", "GLUD1;GLUD2", # Pericentral heps
                  "PCK1", "Saa2;Saa1-1", "HAMP") # Periportal heps
# Alternative: get desired genes for PBMCs

# TAKE ONE involving building dendrogram and Seurat heatmap
sobj <- GetResidual(sobj,
                    features = plotFeatures,
                    assay = "SCT")
# Isolate the scale data values for these genes
featureMatrix <- sobj@assays$SCT@scale.data[plotFeatures,]
# Subsample matrix
subMat <- featureMatrix[,sample(ncol(featureMatrix), size = 10000)]
# Cluster genes in subsampled matrix
dendro <- as.dendrogram(hclust(d = dist(x = subMat))) # Note that the leafs are the row names
# Also try to do dendro for clusters
# First get cluster averages
woodchuckClusterAverages <- AverageExpression(sobj,
                                              assays = "SCT",
                                              slot = "counts")
dendroClusts <- as.dendrogram(hclust(d = dist(x = t(woodchuckClusterAverages$SCT))))
# To visualize dendrogram:
# ggdendro::ggdendrogram(dendroClusts, rotate = TRUE)
Idents(sobj) <- factor(x = Idents(sobj), 
                       levels = colnames(woodchuckClusterAverages$SCT)[unlist(dendroClusts)])
# Reset to increasing order
# Idents(sobj) <- factor(x = Idents(sobj), levels = sort(as.numeric(levels(Idents(sobj)))))
heatmap <- DoHeatmap(subset(sobj, downsample = 500),
                     features = rownames(subMat)[unlist(dendro)])
heatmap
# Try grouping by refined cell labels as an example to see if it works better
heatmap2 <- DoHeatmap(subset(sobj, downsample = 500),
                      features = rownames(subMat)[unlist(dendro)],
                      group.by = "scina_labels_refined",
                      angle = 90, size = 3
                      ) & NoLegend() #& coord_flip()
heatmap2

# HEATMAP TAKE TWO with GGPLOT
# Isolate gene expression matrix
clustAvg <- AverageExpression(sobj,
                              assays = "SCT",
                              slot = "counts",
                              group.by = "scina_labels_refined")
geneMat <- t(scale(t(clustAvg[[1]][plotFeatures,]))) # Scale across genes
cormat <- reshape2::melt(as.matrix(geneMat))
ggplot2::ggplot(cormat, ggplot2::aes(x = Var2, y = Var1, fill = value)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(low = "blue", high = "darkred",
                                name = "Expression value") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                     vjust = 1,
                                                     hjust = 1))

# HEATMAP TAKE THREE with pheatmap
geneMat <- t(clustAvg[[1]][plotFeatures,]) # No scaling
pheatmap::pheatmap(geneMat,
                   scale = "column")

# Now make dotplot using same clustering of features
dotplot <- DotPlot(subset(sobj, downsample = 500),
                   features = rownames(subMat)[unlist(dendro)]) + 
  RotatedAxis()
dotplot
# Now try cell labels
dotplot2 <- DotPlot(subset(sobj, downsample = 500),
                    features = rownames(subMat)[unlist(dendro)],
                    group.by = "scina_labels_refined") + 
  RotatedAxis() +
  theme(axis.text = element_text(size = 8))
dotplot2

# Give woodchuck cells general labels at low resolution


# Correlation functions/set-up: DO FOR ALL

# Read in ortholog table
geneNameTable <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/collectedOrthofinderPairings.tsv",
                            sep = "\t",
                            header = TRUE)
woodchuckClusterAverages <- AverageExpression(sobj,
                                              assays = "SCT",
                                              slot = "scale.data")
# Scale data
woodchuckClusterAverages$SCT <- na.omit(t(scale(t(woodchuckClusterAverages$SCT))))
# Grab gene names from Seurat object
uniqueHier <- row.names(woodchuckClusterAverages$SCT)
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable to get correct order (notice uniqueHier is on left)
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# Get orthologs from either mouse or human
species <- "human"
if (species == "human") {
  # If human one-to-one ortho has NA, replace with mikado_final_sc2_stringent_noMito_protein column
  # This is to avoid and potential mistakes in recognizing things it shouldn't be recognizing
  newNames$speciesOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$uniqueHier, newNames$humanOneToOne)
} else if (species == "mouse") {
  newNames$speciesOneToOne <- ifelse(is.na(newNames$mouseOneToOne), newNames$uniqueHier, newNames$mouseOneToOne)
} else if (species == "woodchuck") {
  newNames$speciesOneToOne <- newNames$uniqueHier
}
# Grab dataframe
woodchuckClusterAverages <- woodchuckClusterAverages$SCT
# Replace names with one-to-one orthologue of particular species
row.names(woodchuckClusterAverages) <- newNames$speciesOneToOne
# Make sure formatted correctly
woodchuckClusterAverages <- as.data.frame(woodchuckClusterAverages)
# Order by gene name
woodchuckClusterAverages <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]

# Correlations of woodchuck to human (PBMCs)

# Now need to read in 68K PBMC data
humanPBMC <- read.csv("~/Dropbox/Zoe/scf_version/analysis/correlationTests/68K_pbmc_data/68K_enrichedGenes.csv",
                      header = FALSE)
# Get rid of top row
humanPBMC <- humanPBMC[-1,]
# Separate all cell and myloid cell data
allCells <- select(humanPBMC, V1, V2, V3)
myeloid <- select(humanPBMC, V5, V6, V7)
# Rename and get rid of first row
colnames(allCells) <- c("Cluster", "Gene", "Enrichment")
colnames(myeloid) <- c("Cluster", "Gene", "Enrichment")
# Get rid of top row
allCells <- allCells[-1,]
myeloid <- myeloid[-1,]
# I can then filter the rows and bind the data frames back together by gene name
for (j in 1:10) {
  dat <- dplyr::filter(allCells, Cluster == j)
  #dat <- get(paste("allCells", j, sep = ""))
  # Multiply enrichment values by -1 because the signs are backwards???
  dat$Enrichment <- as.numeric(dat$Enrichment) * -1
  dat <- dplyr::select(dat, Gene, Enrichment)
  colnames(dat) <- c("Gene", paste("Enrichment", j, sep = ""))
  if (j == 1) {
    allCellsMatrix <- dat
  }
  else {
    allCellsMatrix <- dplyr::full_join(allCellsMatrix, dat, by = "Gene")
  }
}
# Make row names gene names
rownames(allCellsMatrix) <- allCellsMatrix$Gene
allCellsMatrix <- dplyr::select(allCellsMatrix, -Gene)
# Make column names cell types
colnames(allCellsMatrix) <- c("Activated CD8+", "Naive CD8+", "Memory and Reg T",
                              "Naive CD4+", "NK", "CD8+", "B", "Megakaryocytes",
                              "Monocytes and Dendritic", "B, Dendritic, T")
# Scale across columns
allCellsMatrix <- t(scale(t(allCellsMatrix)))
# Order genes alphabetically by gene name
allCellsMatrix <- allCellsMatrix[order(row.names(allCellsMatrix)),]

# Correlations of woodchuck liver to human 5 liver map

# Find cluster averages of human liver data
load("~/Dropbox/Zoe/scf_version/analysis/correlationTests/HumanLiver.RData")
# Run SCTransform
HumanLiverSeurat <- UpdateSeuratObject(HumanLiverSeurat)
HumanLiverSeurat <- SCTransform(HumanLiverSeurat)
humanClusterAverages <- AverageExpression(HumanLiverSeurat,
                                          assays = "SCT",
                                          slot = "scale.data")
# Replace cluster numbers with names
colnames(humanClusterAverages$SCT) <- c("Hep 1", "Alpha-beta T cells", "Hep 2",
                                        "Inflammatory macs", "Hep 3", "Hep 4",
                                        "Plasma cells", "NK-like cells", "Gamma-delta T cells",
                                        "Non-inflammatory macs", "Periportal LSECs", "Central venous LSECs",
                                        "Portal endothelial cells", "Hep 5", "Hep 6",
                                        "Mature B cells", "Cholangiocytes", "Gamma-delta T cells 2",
                                        "Erythroid cells", "Hepatic stellate cells")
# If only looking at specific clusters
#humanClusterAverages$SCT <- humanClusterAverages$SCT[,c("3","1","15","6","14","5")]
# Otherwise go straight to here:
humanClusterAverages$SCT <- na.omit(t(scale(t(humanClusterAverages$SCT))))
# Grab gene names
humanGenes <- row.names(humanClusterAverages$SCT)
# Now turn into large dataframe
allCellsMatrix <- as.data.frame(humanClusterAverages$SCT)
# Order by row name
allCellsMatrix <- allCellsMatrix[order(row.names(allCellsMatrix)),]

# Prep Aizarani data

# Read in Aizarani dataset
aizarani <- readRDS("~/Dropbox/Zoe/scf_version/analysis/correlationTests/GSE124395_Normalhumanliverdata.RData")
# Read in clusters and label cells
aizaraniClusters <- read.table("~/Dropbox/Zoe/scf_version/analysis/correlationTests/GSE124395_clusterpartition.txt")
# Only keep cells in the cluster object
aizarani <- aizarani[,intersect(colnames(aizarani),row.names(aizaraniClusters))]
# Create Seurat object
aizarani <- CreateSeuratObject(counts = aizarani)
# Run SCTransform
aizarani <- SCTransform(aizarani)
# Add cluster IDs
Idents(aizarani) <- aizaraniClusters$sct.cpart
# Get cluster averages
aizaraniAverages <- AverageExpression(aizarani,
                                      assays = "SCT",
                                      slot = "scale.data")
aizaraniAverages$SCT <- na.omit(t(scale(t(aizaraniAverages$SCT))))
# Grab gene names
aizaraniGenes <- row.names(aizaraniAverages$SCT)
# Now turn into large dataframe
allCellsMatrix <- as.data.frame(aizaraniAverages$SCT)
# Order by row name
allCellsMatrix <- allCellsMatrix[order(row.names(allCellsMatrix)),]

# Correlate woodchuck liver with PBMCs

# Start with liver and read in woodchuck PBMCs again
load("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/integrated_PBMC_cca_kanchor5_scClustViz.RData")
Idents(scSeurat) <- "integrated_snn_res.0.6"
# Find cluster averages
pbmcClusterAverages <- AverageExpression(scSeurat,
                                         assays = "SCT",
                                         slot = "scale.data")
pbmcClusterAverages <- as.data.frame(na.omit(t(scale(t(pbmcClusterAverages$SCT)))))
# Order by row name
allCellsMatrix <- pbmcClusterAverages[order(row.names(pbmcClusterAverages)),]

# Now go back to doing for all datasets

speciesData <- "Human Liver (Aizarani et al.)"
woodchuckData <- "Woodchuck Liver"
# Now find intersecting genes
matches <- intersect(row.names(allCellsMatrix),
                     row.names(woodchuckClusterAverages))
# Look at how many genes matched
length(matches)
# Make new matrices with only matching gene names
toCor <- allCellsMatrix[matches,]
woodchuckAveragesCor <- woodchuckClusterAverages[matches,]
# Do Pearson
pearVal <- cor(toCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = paste("Pearson correlation of", speciesData, "vs", woodchuckData),
        xlab = woodchuckData,
        ylab = speciesData)
        #margins = c(6,11))
#Rowv = NA,
#Colv = NA)
# Do Spearman
spearVal <- cor(toCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = paste("Spearman correlation of", speciesData, "vs", woodchuckData),
        xlab = woodchuckData,
        ylab = speciesData)
        #margins = c(6,11))
#Rowv = NA,
#Colv = NA)
