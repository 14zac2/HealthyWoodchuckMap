library(Seurat)
library(dplyr)
library(ggplot2)

# Load object
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_TLH_rpca_kanchor20_mito50_dropletQC_scClustViz.RData")

# Find top markers
Idents(scSeurat) <- "integrated_snn_res.0.4"
markers <- FindAllMarkers(scSeurat, only.pos = TRUE)
# Get top 5
top.markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(scSeurat, features = top.markers$gene, angle = 0, size = 3) & NoLegend() +
  theme(axis.text.y = element_text(size = 5))

# Now download markers
