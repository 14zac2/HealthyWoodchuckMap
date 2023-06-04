library(Seurat)
library(dplyr)
library(scClustViz)
library(viridis)
library(presto)
library(cluster)
library(viridisLite)
library(shiny)

dropletqc <- "no_dropletQC"

# Get woodchuck metadata
source("~/Dropbox/Zoe/scf_version/analysis/scripts/woodchuckMetadata.R")

# Isolate healthy woodchucks
samples <- woodchuck_info[woodchuck_info$tissue == "Healthy",]

# Further isolate if desired
samples <- samples[grepl("^L", samples$woodchuck),]
# For PBMCs
samples <- samples[samples$woodchuck != "3321",]
samples$sample_type <- "PBMC"
samples$prefix <- gsub("BIOPSY|TLH","PBMC", samples$prefix)

objects <- c()

# In a loop, load all dropletQC filtered objects
for (row in 1:nrow(samples)) {
  load(paste("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/",
             dropletqc, "/", samples$woodchuck[row], "/", samples$woodchuck[row], "_", samples$sample_type[row], 
             "_scClustViz.RData", sep = ""))
  assign(paste("sobj", row, sep = ""), scSeurat)
  # Add this to vector
  objects <- c(objects, get(paste("sobj", row, sep = "")))
}

merged_sobj <- merge(x = sobj1, y = objects[2:nrow(samples)],
                     merge.data = TRUE,
                     add.cell.ids = samples$prefix)

# Make variable features that of the first object
VariableFeatures(merged_sobj[["SCT"]]) <- rownames(sobj1[["SCT"]]@scale.data)

# Normalise without SCTransform
merged_sobj <- merged_sobj %>%
  RunPCA() %>% 
  FindNeighbors(dims = 1:30) %>% 
  RunUMAP(dims = 1:30) #%>%
merged_sobj <- RunTSNE(merged_sobj)

# Prep to find markers
merged_sobj <- PrepSCTFindMarkers(merged_sobj, assay = "SCT", verbose = TRUE)

# Remove all objects at this point
rm(list=ls(pattern="^sobj"))
gc()

# Cluster
for (res in seq(from = 0.2, to = 2.0, by = 0.2)) {
  merged_sobj <- FindClusters(merged_sobj, resolution = res)
}
Idents(merged_sobj) <- "SCT_snn_res.0.2"

saveRDS(merged_sobj,
        file = "~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/LSampsMerged_dropletQCFiltered.RDS")

DimPlot(merged_sobj, label = TRUE)
DimPlot(merged_sobj, label = TRUE, group.by = "scina_labels")
DimPlot(merged_sobj, group.by = "orig.ident")
FeaturePlot(merged_sobj, features = c("nuclear_fraction","percent.mt"))
FeaturePlot(merged_sobj, features = c("CYP2E1","FETUB","PCK1","C8G"))

# OR scClustViz
library(scClustViz)
library(viridis)
library(presto)
library(cluster)
library(viridisLite)
library(shiny)
scSeurat <- merged_sobj
DE_bw_clust <- TRUE
FDRthresh <- 0.01
seurat_resolution <- 0
sCVdata_list <- list()

while(DE_bw_clust) {
  seurat_resolution <- seurat_resolution + 0.2
  # ^ iteratively incrementing resolution parameter
  
  scSeurat <- FindClusters(scSeurat,
                           resolution=seurat_resolution,
                           print.output=F)
  message(" ")
  message("------------------------------------------------------")
  message(paste0("--------  res.",seurat_resolution," with ",
                 length(levels(Idents(scSeurat)))," clusters --------"))
  message("------------------------------------------------------")
  curr_sCVdata <- CalcSCV(inD=scSeurat,
                          cl=Idents(scSeurat),
                          # ^ your most recent clustering results get stored in the Seurat "ident" slot
                          assayType = "SCT",
                          exponent=exp(1),
                          pseudocount=1,
                          DRthresh=0.1,
                          DRforClust="pca",
                          calcSil=T,
                          calcDEvsRest=T,
                          calcDEcombn=T)
  
  DE_bw_NN <- sapply(DEneighb(curr_sCVdata,FDRthresh),nrow)
  # ^ counts # of DE genes between neighbouring clusters at your selected FDR threshold
  
  if (min(DE_bw_NN) < 1 | seurat_resolution > 2) { DE_bw_clust <- FALSE }
  # ^ If no DE genes between nearest neighbours, don't loop again.
  
  sCVdata_list[[paste0("res.",seurat_resolution)]] <- curr_sCVdata
  gc()
}

# Save scClustViz object
save(sCVdata_list, scSeurat,
     file = "~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/merged_PBMC_scClustViz.RData")
