# Script to find markers intersecting between datasets

# Single-cell dataset
sc <- readRDS("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/allIntegrated_cca_kanchor5_noBiopsyHeps_dropletQCFiltered.RDS")
Idents(sc) <- "integrated_snn_res.1.4"

# PBMC dataset
load("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/integrated_PBMC_cca_kanchor5_scClustViz.RData")
sc <- scSeurat
Idents(sc) <- "integrated_snn_res.0.6"

# Single-nuc dataset
#load("~/Dropbox/Zoe/scf_version/analysis/multiome/healthy_integration/L192_L212_multiome_onlyRNA_rpca_dropletQC_0.5nf_scClustViz.RData")
#sn <- scSeurat
#Idents(sn) <- "integrated_snn_res.0.6"

# Look at markers for an individual dataset that pass a certain FDR threshold - SEURAT method
sc <- readRDS("~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/allIntegrated_cca_kanchor5_noBiopsyHeps_dropletQCFiltered.RDS")
Idents(sc) <- "integrated_snn_res.1.4"

markers <- FindAllMarkers(sc,
                          slot = "scale.data",
                          assay = "SCT",
                          only.pos = TRUE,
                          recorrect_umi = FALSE)
sigMarkers <- markers[markers$p_val_adj < 0.05,]
# Check how many markers are in dataset
for (clust in levels(as.factor(sigMarkers$cluster))) {
  print(paste("Cluster", clust, "has", sum(sigMarkers$cluster == clust), "DE genes."))
}
# Print top genes for each cluster
for (clust in levels(as.factor(sigMarkers$cluster))) {
  cat("Cluster ", clust, ": ", sep = "")
  cat(sigMarkers$gene[sigMarkers$cluster == clust][1:25], "\n", sep = ", ")
}
write.table(sigMarkers,
            #file = "~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/dropletQC_filtered/allMarkers_allIntegrated_cca_kanchor5_noBiopsyHeps_dropletQCFiltered_res1.4.tsv",
            file = "~/Dropbox/Zoe/scf_version/analysis/healthy_sc/seurat_objects/no_dropletQC/allMarkers_integrated_PBMC_cca_kanchor5_res0.6.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# Choose clusters to compare
sc_cluster <- c(12)
sn_cluster <- c(8)

# Find single-cell markers
sc_genes <- FindMarkers(sc, ident.1 = sc_cluster,
            slot = "scale.data", assay = "SCT",
            only.pos = TRUE, recorrect_umi = FALSE)
paste("There are", nrow(sc_genes), "marker genes.")
sc_genes <- sc_genes[grep("mikado", row.names(sc_genes), invert = TRUE),]
cat(row.names(sc_genes), sep = ", ")

# Find single-nuc markers
#sn_genes <- FindMarkers(sn, ident.1 = sn_cluster,
#                        slot = "scale.data", assay = "SCT",
#                        only.pos = TRUE, recorrect_umi = FALSE)
#paste("There are", nrow(sn_genes), "marker genes.")
#sn_genes <- sn_genes[grep("mikado", row.names(sn_genes), invert = TRUE),]
#cat(row.names(sn_genes), sep = ", ")

# Find those that intersect
#sc_sn_genes <- intersect(rownames(sc_genes), rownames(sn_genes))
#paste("There are", length(sc_sn_genes), "intersecting genes.")
#cat(sc_sn_genes, sep = ", ")

## Now compare to human

# Read in ortholog table
geneNameTable <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/collectedOrthofinderPairings.tsv",
                            sep = "\t",
                            header = TRUE)

# Compare to human single-cell
load("~/Dropbox/Zoe/scf_version/analysis/correlationTests/HumanLiver.RData")
# Run SCTransform
HumanLiverSeurat <- SCTransform(HumanLiverSeurat)
# Replace cluster numbers with names
levels(Idents(HumanLiverSeurat)) <- c("Hep 1", "Alpha-beta T cells", "Hep 2",
                                      "Inflammatory macs", "Hep 3", "Hep 4",
                                      "Plasma cells", "NK-like cells", "Gamma-delta T cells",
                                      "Non-inflammatory macs", "Periportal LSECs", "Central venous LSECs",
                                      "Portal endothelial cells", "Hep 5", "Hep 6",
                                      "Mature B cells", "Cholangiocytes", "Gamma-delta T cells 2",
                                      "Erythroid cells", "Hepatic stellate cells")
# Choose clusters to compare
woodchuck_cluster <- c(8)
human_cluster <- c("Plasma cells")
# Find woodchuck single-cell markers
woodchuck_genes <- FindMarkers(sc, ident.1 = woodchuck_cluster,
                        slot = "scale.data", assay = "SCT",
                        only.pos = TRUE, recorrect_umi = FALSE)
paste("There are", nrow(woodchuck_genes), "marker genes.")
woodchuck_genes <- woodchuck_genes[grep("mikado", row.names(woodchuck_genes), invert = TRUE),]
# Find human single-cell markers
human_genes <- FindMarkers(HumanLiverSeurat, ident.1 = human_cluster,
                               slot = "scale.data", assay = "SCT",
                               only.pos = TRUE, recorrect_umi = FALSE)
paste("There are", nrow(human_genes), "marker genes.")
# Grab woodchuck gene names
uniqueHier <- row.names(woodchuck_genes)
# Replace these gene names with the human name
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# If mikado protein column is blank, replace with uniqueHier (non-coding genes)
newNames$mikado_final_sc2_stringent_noMito_protein <- ifelse(is.na(newNames$mikado_final_sc2_stringent_noMito_protein), newNames$uniqueHier, newNames$mikado_final_sc2_stringent_noMito_protein)
# If human has NA, replace with mikado name column
newNames$humanOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$humanOneToOne)
# Add gene names as human
row.names(woodchuck_genes) <- newNames$humanOneToOne
# Find those that intersect
woodchuck_human_genes <- intersect(rownames(woodchuck_genes), rownames(human_genes))
paste("There are", length(woodchuck_human_genes), "intersecting genes.")
cat(woodchuck_human_genes, sep = ", ")

# Compare to human single-nuc
humanNuc <- readRDS("~/Dropbox/Zoe/scf_version/analysis/correlationTests/single_nuc_20_human_map.rds")
# Isolate single-nuc only
humanNuc <- subset(humanNuc, subset = assay_type == "single_nuc")
# Run SCTransform
#options(future.globals.maxSize = 8000 * 1024^2)
humanNuc <- SCTransform(humanNuc)
# Change idents
Idents(humanNuc) <- humanNuc@meta.data$Manual_Annotation
# See human single-nuc cell types
levels(Idents(humanNuc))
# Choose clusters to compare
woodchuck_sn_clust <- c(8)
human_sn_clust <- c("Bcells")
# Find woodchuck single-cell markers
woodchuck_sn_genes <- FindMarkers(sn, ident.1 = woodchuck_sn_clust,
                               slot = "scale.data", assay = "SCT",
                               only.pos = TRUE, recorrect_umi=FALSE)
paste("There are", nrow(woodchuck_sn_genes), "marker genes.")
woodchuck_sn_genes <- woodchuck_sn_genes[grep("mikado", row.names(woodchuck_sn_genes), invert = TRUE),]
# Find human single-cell markers
human_sn_genes <- FindMarkers(humanNuc, ident.1 = human_sn_clust,
                           slot = "scale.data", assay = "SCT",
                           only.pos = TRUE, recorrect_umi=FALSE)
paste("There are", nrow(human_sn_genes), "marker genes.")
# Grab woodchuck gene names
uniqueHier <- row.names(woodchuck_sn_genes)
# Replace these gene names with the human name
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# If mikado protein column is blank, replace with uniqueHier (non-coding genes)
newNames$mikado_final_sc2_stringent_noMito_protein <- ifelse(is.na(newNames$mikado_final_sc2_stringent_noMito_protein), newNames$uniqueHier, newNames$mikado_final_sc2_stringent_noMito_protein)
# If human has NA, replace with mikado name column
newNames$humanOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$humanOneToOne)
# Add gene names as human
row.names(woodchuck_sn_genes) <- newNames$humanOneToOne
# Find those that intersect
woodchuck_human_sn_genes <- intersect(rownames(woodchuck_sn_genes), rownames(human_sn_genes))
paste("There are", length(woodchuck_human_sn_genes), "intersecting genes.")
cat(woodchuck_human_sn_genes, sep = ", ")

# Find genes that intersect with all four datasets
# Just look for interactions between human-woodchuck comparisons
all_intersecting_genes <- intersect(woodchuck_human_genes, woodchuck_human_sn_genes)
paste("There are", length(all_intersecting_genes), "intersecting genes.")
cat(all_intersecting_genes, sep = ", ")
