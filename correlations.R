# Pearson and Spearman correlations against human and mouse liver
library(scClustViz)
library(Seurat)
library(dplyr)
library(Hmisc)

# Load scClustViz object
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_TLH_rpca_kanchor20_mito50_dropletQC_scClustViz.RData")
#load("~/Dropbox/Zoe/scf_version/analysis/multiome/healthy_integration/L192_L212_multiome_onlyRNA_rpca_dropletQC_0.5nf_scClustViz.RData")

# Read in ortholog table
geneNameTable <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/collectedOrthofinderPairings.tsv",
                            sep = "\t",
                            header = TRUE)

# Find cluster averages of woodchuck data
# Use SCTransform output
# Choose clustering resolution
# Idents(scSeurat) <- "SCT_snn_res.0.8"
Idents(scSeurat) <- "integrated_snn_res.0.4"
# Looks like there are counts but not scaled data for every gene
woodchuckClusterAverages <- AverageExpression(scSeurat,
                                              assays = "SCT",
                                              slot = "scale.data")
# If only looking at certain clusters
#woodchuckClusterAverages$SCT <- woodchuckClusterAverages$SCT[,c("0","1","2","4","5","6","9","10","11")]
# Otherwise go straight to here:
woodchuckClusterAverages$SCT <- na.omit(t(scale(t(woodchuckClusterAverages$SCT))))
# Grab gene names
uniqueHier <- row.names(woodchuckClusterAverages$SCT)
# Replace these gene names with the human name
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# If mikado protein column is blank, replace with uniqueHier (non-coding genes)
newNames$mikado_final_sc2_stringent_noMito_protein <- ifelse(is.na(newNames$mikado_final_sc2_stringent_noMito_protein), newNames$uniqueHier, newNames$mikado_final_sc2_stringent_noMito_protein)
# If human or mouse has NA, replace with mikado name column
newNames$humanOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$humanOneToOne)
newNames$mouseOneToOne <- ifelse(is.na(newNames$mouseOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$mouseOneToOne)
# Fix format
woodchuckClusterAverages <- woodchuckClusterAverages$SCT
# Add gene names as human
row.names(woodchuckClusterAverages) <- newNames$humanOneToOne
# Make sure formatted as list of dataframes
woodchuckClusterAverages <- as.data.frame(woodchuckClusterAverages)

# Find cluster averages of human liver data
load("~/Dropbox/Zoe/scf_version/analysis/correlationTests/HumanLiver.RData")
# Run SCTransform
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
humanClusterAverages <- as.data.frame(humanClusterAverages$SCT)
# Order by row name
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
humanAveragesOrdered <- humanClusterAverages[order(row.names(humanClusterAverages)),]
# Now find intersecting genes
humanMatches <- intersect(row.names(woodchuckAveragesOrdered),
                          row.names(humanAveragesOrdered))
# Look at how many genes matched
length(humanMatches)
# Make new matrices with only matching gene names
humanAveragesCor <- humanAveragesOrdered[humanMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[humanMatches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(humanAveragesCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(humanAveragesCor)[, i], as.matrix(woodchuckAveragesCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(humanAveragesCor)),
       dimnames=list(colnames(as.matrix(humanAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(humanAveragesCor)),
       dimnames=list(colnames(as.matrix(humanAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Pearson correlation of woodchuck vs human",
        xlab = "Woodchuck single-cell clusters",
        ylab = "Human single-cell clusters",
        margins = c(5,11))
        #cexRow = 3, cexCol = 3,
        #Rowv = NA)
pearVal <- cor(humanAveragesCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of woodchuck vs human",
        xlab = "Woodchuck TLH clusters",
        ylab = "Human liver clusters",
        #cexRow = 3, cexCol = 3,
        Rowv = NA)
# Do Spearman
spearVal <- cor(humanAveragesCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of woodchuck vs human (MacParland et al. 2018)",
        xlab = "Woodchuck TLH clusters")
       # ylab = "Human liver clusters")
        #cexRow = 3, cexCol = 3,
        #Rowv = NA)
# Do human against itself
pearVal <- cor(humanAveragesCor, humanAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of human vs human",
        xlab = "Human liver clusters",
        ylab = "Human liver clusters")

# Look at mouse values
mouseLiverTest <- read.csv("~/Dropbox/Zoe/scf_version/analysis/correlationTests/mouseLiverTests.csv",
                           row.names = 1,
                           header = TRUE)
# Count number of genes
nrow(mouseLiverTest)
mouseLiverTest <- na.omit(mouseLiverTest)
nrow(mouseLiverTest)
# Note that the second zone values are error or something (s.e.m.)
# Isolate p vals less than 1 x 10x-25 from Kruskal Wallis
mouseLiverTest <- dplyr::filter(mouseLiverTest, p.values < 1*10^-02) # this was -25
# Or q vals
#mouseLiverTest <- dplyr::filter(mouseLiverTest, q.values < 0.01)
# Check number of genes remaining
nrow(mouseLiverTest)
# Now only keep desired columns
mouseClusterAverages <- dplyr::select(mouseLiverTest, Layer.1:Layer.9)
# Scale across columns
mouseClusterAverages <- t(scale(t(mouseClusterAverages)))
# Genes are already ordered alphabetically by gene name
mouseAveragesOrdered <- mouseClusterAverages[order(row.names(mouseClusterAverages)),]
# Get woodchuck genes for mouse comparison
row.names(woodchuckClusterAverages) <- newNames$mouseOneToOne
# Reorder for mouse
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]

# Now find intersecting genes
mouseMatches <- intersect(row.names(woodchuckAveragesOrdered),
                          row.names(mouseAveragesOrdered))
# Look at how many genes matched
length(mouseMatches)
# Make new matrices with only matching gene names
mouseAveragesCor <- mouseAveragesOrdered[mouseMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[mouseMatches,]
colnames(mouseAveragesCor) <- c("1","2","3","4","5","6","7","8","9")
# Do spearman with p-val
tmp <- with(expand.grid(seq(ncol(as.matrix(mouseAveragesCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(mouseAveragesCor)[, i], as.matrix(woodchuckAveragesCor)[, j],
                                           method = "spearman"),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(mouseAveragesCor)),
               dimnames=list(colnames(as.matrix(mouseAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(mouseAveragesCor)),
                dimnames=list(colnames(as.matrix(mouseAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Spearman of woodchuck vs mouse",
        xlab = "Woodchuck TLH clusters",
        ylab = "Mouse liver zones",
        Rowv = NA)
# Do Pearson
pearVal <- cor(mouseAveragesCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of woodchuck vs mouse",
        xlab = "Woodchuck TLH clusters",
        ylab = "Mouse liver zones",
        Rowv = NA)
        #Colv = NA,
        #cexRow = 3,
        #cexCol = 3)
# Do Spearman
spearVal <- cor(mouseAveragesCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of woodchuck vs mouse",
        xlab = "Woodchuck TLH clusters",
        ylab = "Mouse liver zones",
        Rowv = NA)
        #Colv = NA,
        #cexCol = 2,
        #cexRow = 3)

# Can also try other mouse one
# Periportal
mouseZone1 <- read.csv("~/Dropbox/Zoe/woodchuck_si_aug18/analysis/correlationTests/Mouse3ZonePaper/Zone1hepVSrest.csv",
                           row.names = 1,
                           header = TRUE)
# Midzonal
mouseZone2 <- read.csv("~/Dropbox/Zoe/woodchuck_si_aug18/analysis/correlationTests/Mouse3ZonePaper/Zone2hepVSrest.csv",
                       row.names = 1,
                       header = TRUE)
# Pericentral
mouseZone3 <- read.csv("~/Dropbox/Zoe/woodchuck_si_aug18/analysis/correlationTests/Mouse3ZonePaper/Zone3hepVSrest.csv",
                       row.names = 1,
                       header = TRUE)
# For each, filter for low adjusted p-vals and then remove that row
for (j in 1:3) {
        dat <- get(paste("mouseZone", j, sep = ""))
        dat <- dplyr::filter(dat, adj_pval < 0.05)
        dat <- dplyr::select(dat, logFC)
        colnames(dat) <- paste("logFC", j, sep = "_")
        dat$gene <- row.names(dat)
        if (j == 1) {
                allMouseZones <- dat
        }
        #assign(paste("mouseZone", j, sep = ""), dat)
        else {
                allMouseZones <- inner_join(allMouseZones, dat, by = "gene")
        }
}
# 1573 in common
# Rename with rows as genes
row.names(allMouseZones) <- allMouseZones$gene
# Now get rid of gene row
allMouseZones <- dplyr::select(allMouseZones, -gene)
# Scale across columns
allMouseZones <- t(scale(t(allMouseZones)))
# Order by gene name
allMouseZones <- allMouseZones[order(row.names(allMouseZones)),]
# Get woodchuck genes for mouse comparison
row.names(woodchuckClusterAverages) <- newNames$mouseOneToOne
# Reorder for mouse
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
# Now find intersecting genes
allZonesMatches <- intersect(row.names(woodchuckAveragesOrdered),
                             row.names(allMouseZones))
# Look at how many genes matched
length(allZonesMatches) # 1149
# Make new matrices with only matching gene names
allZonesCor <- allMouseZones[allZonesMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[allZonesMatches,]
colnames(allZonesCor) <- c("Periportal","Midzonal","Pericentral")
# Reorder so matches 7 mouse zone order
allZonesCor <- dplyr::select(as.data.frame(allZonesCor), Pericentral, Midzonal, Periportal)
# Do Pearson
pearVal <- cor(allZonesCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of woodchuck vs 3 mouse zones",
        xlab = "Woodchuck TLH clusters",
        ylab = "Mouse liver zones",
        Rowv = NA)
        #Colv = NA)
# Do Spearman
spearVal <- cor(allZonesCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of woodchuck vs 3 mouse zones",
        xlab = "Woodchuck TLH clusters",
        ylab = "Mouse liver zones",
        Rowv = NA)
        #Colv = NA

# Try human
humanMouseOrder <- humanClusterAverages[order(row.names(humanClusterAverages)),]
# Make lower for mouse match
row.names(humanMouseOrder) <- tolower(row.names(humanMouseOrder))
# Make lower case in mouse, too
mouseLower <- mouseAveragesOrdered
row.names(mouseLower) <- tolower(row.names(mouseLower))
# Find matches
humanMouseMatches <- intersect(row.names(humanMouseOrder),
                          row.names(mouseLower))
length(humanMouseMatches)
# Make new matrices with only matching gene names
mouseHumAveragesCor <- mouseLower[humanMouseMatches,]
humMouseAveragesCor <- humanMouseOrder[humanMouseMatches,]
# Do Pearson
pearVal <- cor(mouseHumAveragesCor, humMouseAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of mouse vs human",
        xlab = "Mouse liver zones",
        ylab = "Human liver clusters",
        Rowv = NA,
        Colv = NA)
# Do Spearman
spearVal <- cor(mouseHumAveragesCor, humMouseAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of mouse vs human",
        xlab = "Mouse liver zones",
        ylab = "Human liver clusters",
        Rowv = NA,
        Colv = NA)

# Now try other mouse zone paper


# Also trying for PBMCs
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_PBMC_rpca_scClustViz.RData")
Idents(scSeurat) <- "integrated_snn_res.0.6"
# Get average gene expression across clusters
woodchuckClusterAverages <- AverageExpression(scSeurat,
                                              assays = "SCT",
                                              slot = "scale.data")
# Scale data
woodchuckClusterAverages$SCT <- na.omit(t(scale(t(woodchuckClusterAverages$SCT))))
# Grab gene names
uniqueHier <- row.names(woodchuckClusterAverages$SCT)
# Replace these gene names with the human ortholog
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable to get correct order (notice uniqueHier is on left)
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# If human one-to-one ortho has NA, replace with mikado_final_sc2_stringent_noMito_protein column
# This is to avoid and potential mistakes in recognizing things it shouldn't be recognizing
newNames$humanOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$uniqueHier, newNames$humanOneToOne)
# Grab dataframe
woodchuckClusterAverages <- woodchuckClusterAverages$SCT
# Replace names
row.names(woodchuckClusterAverages) <- newNames$humanOneToOne
# Make sure formatted correctly
woodchuckClusterAverages <- as.data.frame(woodchuckClusterAverages)

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
woodchuckClusterAverages <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
# Now find intersecting genes
pbmcMatches <- intersect(row.names(allCellsMatrix),
                                   row.names(woodchuckClusterAverages))
# Look at how many genes matched
length(pbmcMatches) # 477
# Make new matrices with only matching gene names
humanPbmcCor <- allCellsMatrix[pbmcMatches,]
woodchuckAveragesCor <- woodchuckClusterAverages[pbmcMatches,]
# Do Pearson
pearVal <- cor(humanPbmcCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of human 68K PBMCs vs woodchuck PBMCs",
        xlab = "Woodchuck clusters",
        ylab = "Human clusters",
        margins = c(6,11))
        #Rowv = NA,
        #Colv = NA)
# Do Spearman
spearVal <- cor(humanPbmcCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of human vs woodchuck PBMCs",
        xlab = "Woodchuck clusters",
        ylab = "Human clusters")
        #Rowv = NA,
        #Colv = NA)

# Compare woodchucks to eachother
load("L202/TLH/L202_TLH_cellranger_mito50_scClustViz.RData")
L202 <- scSeurat
load("L212/TLH/L212_TLH_cellranger_mito50_scClustViz.RData")
L212 <- scSeurat
# Set clustering resolutions
Idents(L202) <- "SCT_snn_res.0.2"
Idents(L212) <- "SCT_snn_res.0.8"
# Get average cluster expression
L202Av <- AverageExpression(L202, assays = "SCT", slot = "scale.data")
L212Av <- AverageExpression(L212, assays = "SCT", slot = "scale.data")
L202Av <- as.data.frame(na.omit(t(scale(t(L202Av$SCT)))))
L212Av <- as.data.frame(na.omit(t(scale(t(L212Av$SCT)))))
# Order by row name
L202AvOrd <- L202Av[order(row.names(L202Av)),]
L212AvOrd <- L212Av[order(row.names(L212Av)),]
# Now find intersecting genes
woodMatches <- intersect(row.names(L202AvOrd), row.names(L212AvOrd))
# Look at how many genes matched
length(woodMatches)
# Make new matrices with only matching gene names
L202AvCor <- L202AvOrd[woodMatches,]
L212AvCor <- L212AvOrd[woodMatches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(L202AvCor))), seq(ncol(as.matrix(L212AvCor)))),
            mapply(function(i, j) cor.test(as.matrix(L202AvCor)[, i], as.matrix(L212AvCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(L202AvCor)),
               dimnames=list(colnames(as.matrix(L202AvCor)), colnames(as.matrix(L212AvCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(L202AvCor)),
                dimnames=list(colnames(as.matrix(L202AvCor)), colnames(as.matrix(L212AvCor))))
heatmap(cors,
        main = "Pearson of L202 vs L212",
        xlab = "L212 TLH clusters",
        ylab = "L202 TLH clusters")

# Hibernation correlations
hibMat <- read.csv("~/Dropbox/Zoe/scf_version/analysis/correlationTests/livHib_PRJNA138223.csv")
# Only keep columns of interest
hibMat <- dplyr::select(hibMat, ORF, GSM688560:GSM688571)
# Relabel to make more sense
colnames(hibMat) <- c("gene",
                      "activeLiver1", "activeLiver2", "activeLiver3",
                      "activeLiver4", "activeLiver5", "activeLiver6",
                      "torporLiver1", "torporLiver2", "torporLiver3",
                      "torporLiver4", "torporLiver5", "torporLiver6")
# Calculate averages, if desired
hibMat <- dplyr::mutate(hibMat, active = rowMeans(cbind(activeLiver1, activeLiver2, activeLiver3,
                                                        activeLiver4, activeLiver5, activeLiver6), na.rm=T))
hibMat <- dplyr::mutate(hibMat, torpor = rowMeans(cbind(torporLiver1, torporLiver2, torporLiver3,
                                                        torporLiver4, torporLiver5, torporLiver6), na.rm=T))
hibMat <- dplyr::select(hibMat, gene, active, torpor)
# Have to unique by gene name
dim(hibMat)
hibMat <- dplyr::distinct(hibMat, gene, .keep_all = TRUE)
dim(hibMat)
# Make gene names as row names
row.names(hibMat) <- hibMat$gene
# Now remove gene column
hibMat <- dplyr::select(hibMat, -gene)
# Now turn into large dataframe
hibMat <- as.data.frame(na.omit(t(scale(t(hibMat)))))
#hibMat <- as.data.frame(hibMat)
# Do same processing to woodchuck as would to human comparisons
# Order by row name
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
hibMatOrdered <- hibMat[order(row.names(hibMat)),]
# Now find intersecting genes
hibMatches <- intersect(row.names(woodchuckAveragesOrdered),
                          row.names(hibMatOrdered))
# Look at how many genes matched
length(hibMatches)
# Make new matrices with only matching gene names
hibMatCor <- hibMatOrdered[hibMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[hibMatches,]
# Do correlation test
tmp <- with(expand.grid(seq(ncol(as.matrix(hibMatCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(hibMatCor)[, i], as.matrix(woodchuckAveragesCor)[, j], method = "spearman"),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(hibMatCor)),
               dimnames=list(colnames(as.matrix(hibMatCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(hibMatCor)),
                dimnames=list(colnames(as.matrix(hibMatCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Spearman of woodchuck vs hibernating bear",
        xlab = "Woodchuck TLH clusters",
        ylab = "Hibernating bear")

# Second hibernation (like first)
# Hibernation correlations
hibMat <- read.csv("~/Dropbox/Zoe/scf_version/analysis/correlationTests/livHib_PRJNA114885.csv")
# Only keep columns of interest
hibMat <- dplyr::select(hibMat, Gene_symbol, liver_hibernating_N03.056:liver_active_N05.132)
# Calculate averages, if desired
#hibMat <- dplyr::mutate(hibMat, active = rowMeans(cbind(activeLiver1, activeLiver2, activeLiver3,
#                                                        activeLiver4, activeLiver5, activeLiver6), na.rm=T))
#hibMat <- dplyr::mutate(hibMat, torpor = rowMeans(cbind(torporLiver1, torporLiver2, torporLiver3,
#                                                        torporLiver4, torporLiver5, torporLiver6), na.rm=T))
#hibMat <- dplyr::select(hibMat, gene, active, torpor)
# Have to unique by gene name
dim(hibMat)
hibMat <- dplyr::distinct(hibMat, Gene_symbol, .keep_all = TRUE)
dim(hibMat)
# Make gene names as row names
row.names(hibMat) <- hibMat$Gene_symbol
# Now remove gene column
hibMat <- dplyr::select(hibMat, -Gene_symbol)
# Now turn into large dataframe
hibMat <- as.data.frame(na.omit(t(scale(t(hibMat)))))
#hibMat <- as.data.frame(hibMat)
# Do same processing to woodchuck as would to human comparisons
# Order by row name
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
hibMatOrdered <- hibMat[order(row.names(hibMat)),]
# Now find intersecting genes
hibMatches <- intersect(row.names(woodchuckAveragesOrdered),
                        row.names(hibMatOrdered))
# Look at how many genes matched
length(hibMatches)
# Make new matrices with only matching gene names
hibMatCor <- hibMatOrdered[hibMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[hibMatches,]
# Do correlation test
tmp <- with(expand.grid(seq(ncol(as.matrix(hibMatCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(hibMatCor)[, i], as.matrix(woodchuckAveragesCor)[, j], method = "spearman"),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(hibMatCor)),
               dimnames=list(colnames(as.matrix(hibMatCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(hibMatCor)),
                dimnames=list(colnames(as.matrix(hibMatCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Spearman of woodchuck vs hibernating bear",
        xlab = "Woodchuck TLH clusters",
        ylab = "Hibernating bear",
        Rowv = FALSE)

# More hibernation but also with arousal
hibMat <- read.csv("~/Dropbox/Zoe/scf_version/analysis/correlationTests/livHib_PRJNA168447.csv")
# Only keep columns of interest
hibMat <- dplyr::select(hibMat, -ID_REF, -EST_ID)
# Relabel to make more sense
colnames(hibMat) <- c("gene",
                      "earlyArousal1", "earlyArousal2", "earlyArousal3", "earlyArousal4",
                      "lateArousal1", "lateArousal2", "lateArousal3", "lateArousal4",
                      "earlyTorpor1", "earlyTorpor2", "earlyTorpor3", "earlyTorpor4",
                      "lateTorpor1", "lateTorpor2", "lateTorpor3", "lateTorpor4", "lateTorpor5",
                      "summerActive1", "summerActive2", "summerActive3", "summerActive4", "summerActive5", "summerActive6", "summerActive7")
# Calculate averages, if desired
hibMat <- dplyr::mutate(hibMat, earlyArousal = rowMeans(cbind(earlyArousal1, earlyArousal2, earlyArousal3, earlyArousal4), na.rm=T))
hibMat <- dplyr::mutate(hibMat, lateArousal = rowMeans(cbind(lateArousal1, lateArousal2, lateArousal3, lateArousal4), na.rm=T))
hibMat <- dplyr::mutate(hibMat, earlyTorpor = rowMeans(cbind(earlyTorpor1, earlyTorpor2, earlyTorpor3, earlyTorpor4), na.rm=T))
hibMat <- dplyr::mutate(hibMat, lateTorpor = rowMeans(cbind(lateTorpor1, lateTorpor2, lateTorpor3, lateTorpor4,
                                                            lateTorpor5), na.rm=T))
hibMat <- dplyr::mutate(hibMat, summerActive = rowMeans(cbind(summerActive1, summerActive2, summerActive3, summerActive4,
                                                              summerActive5, summerActive6, summerActive7), na.rm=T))
hibMat <- dplyr::select(hibMat, gene, earlyArousal, lateArousal,
                        earlyTorpor, lateTorpor, summerActive)
# Have to unique by gene name
dim(hibMat)
hibMat <- dplyr::distinct(hibMat, gene, .keep_all = TRUE)
dim(hibMat)
# Make gene names as row names
row.names(hibMat) <- hibMat$gene
# Now remove gene column
hibMat <- dplyr::select(hibMat, -gene)
# Now turn into large dataframe
hibMat <- as.data.frame(na.omit(t(scale(t(hibMat)))))
#hibMat <- as.data.frame(hibMat)
# Do same processing to woodchuck as would to human comparisons
# Order by row name
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
hibMatOrdered <- hibMat[order(row.names(hibMat)),]
# Now find intersecting genes
hibMatches <- intersect(row.names(woodchuckAveragesOrdered),
                        row.names(hibMatOrdered))
# Look at how many genes matched
length(hibMatches)
# Make new matrices with only matching gene names
hibMatCor <- hibMatOrdered[hibMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[hibMatches,]
# Do correlation test
tmp <- with(expand.grid(seq(ncol(as.matrix(hibMatCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(hibMatCor)[, i], as.matrix(woodchuckAveragesCor)[, j], method = "spearman"),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(hibMatCor)),
               dimnames=list(colnames(as.matrix(hibMatCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(hibMatCor)),
                dimnames=list(colnames(as.matrix(hibMatCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Spearman of woodchuck vs hibernating bear",
        xlab = "Woodchuck TLH clusters",
        ylab = "Hibernating bear")

# Single-nuc correlations
# Read in woodchuck data
load("~/Dropbox/Zoe/scf_version/analysis/multiome/healthy_integration/L192_L212_multiome_onlyRNA_rpca_dropletQC_0.5nf_scClustViz.RData")
# Read in ortholog table
geneNameTable <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/collectedOrthofinderPairings.tsv",
                            sep = "\t",
                            header = TRUE)
# Choose clustering resolution
Idents(scSeurat) <- "integrated_snn_res.0.6"
# Get averages from variable genes
woodchuckClusterAverages <- AverageExpression(scSeurat,
                                              assays = "SCT",
                                              slot = "scale.data")
# Use all clusters
woodchuckClusterAverages$SCT <- na.omit(t(scale(t(woodchuckClusterAverages$SCT))))
# Grab gene names
uniqueHier <- row.names(woodchuckClusterAverages$SCT)
# Replace these gene names with the human name
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# If mikado protein column is blank, replace with uniqueHier (non-coding genes)
newNames$mikado_final_sc2_stringent_noMito_protein <- ifelse(is.na(newNames$mikado_final_sc2_stringent_noMito_protein), newNames$uniqueHier, newNames$mikado_final_sc2_stringent_noMito_protein)
# If human or mouse has NA, replace with mikado name column
newNames$humanOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$humanOneToOne)
newNames$mouseOneToOne <- ifelse(is.na(newNames$mouseOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$mouseOneToOne)
# Fix format
woodchuckClusterAverages <- woodchuckClusterAverages$SCT
# Add gene names as human
row.names(woodchuckClusterAverages) <- newNames$humanOneToOne
# Make sure formatted as list of dataframes
woodchuckClusterAverages <- as.data.frame(woodchuckClusterAverages)
# Now read in single-nuc human data
humanNuc <- readRDS("~/Dropbox/Zoe/scf_version/analysis/correlationTests/single_nuc_20_human_map.rds")
# Isolate single-nuc only
humanNuc <- subset(humanNuc, subset = assay_type == "single_nuc")
# Run SCTransform
humanNuc <- SCTransform(humanNuc)
# Change idents
Idents(humanNuc) <- humanNuc@meta.data$sub_annotation
Idents(humanNuc) <- humanNuc@meta.data$Manual_Annotation
# Get cluster averages
humanNucClustAverages <- AverageExpression(humanNuc,
                                           assays = "SCT",
                                           slot = "scale.data")
humanNucClustAverages$SCT <- na.omit(t(scale(t(humanNucClustAverages$SCT))))
# Grab gene names
nucGenes <- row.names(humanNucClustAverages$SCT)
# Now turn into large dataframe
humanNucClustAverages <- as.data.frame(humanNucClustAverages$SCT)
# Order by row name
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
humanAveragesOrdered <- humanNucClustAverages[order(row.names(humanNucClustAverages)),]
# Now find intersecting genes
humanMatches <- intersect(row.names(woodchuckAveragesOrdered),
                          row.names(humanAveragesOrdered))
# Look at how many genes matched
length(humanMatches)
# Make new matrices with only matching gene names
humanAveragesCor <- humanAveragesOrdered[humanMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[humanMatches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(humanAveragesCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(humanAveragesCor)[, i], as.matrix(woodchuckAveragesCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(humanAveragesCor)),
               dimnames=list(colnames(as.matrix(humanAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(humanAveragesCor)),
                dimnames=list(colnames(as.matrix(humanAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Pearson correlation of woodchuck single-nuc vs human single-nuc",
        xlab = "Woodchuck single-nuc clusters",
        ylab = "Human liver clusters",
        margins = c(5,12))
        #cexRow = 3, cexCol = 3,
        #Rowv = NA)
pearVal <- cor(humanAveragesCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of woodchuck single-nuc vs human single-nuc",
        xlab = "Woodchuck TLH clusters",
        ylab = "Human liver clusters",
        margins = c(5,12))
        #cexRow = 3, cexCol = 3,
        #Rowv = NA)
# Do Spearman
spearVal <- cor(humanAveragesCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of woodchuck single-nuc vs human single-nuc",
        xlab = "Woodchuck TLH clusters",
        ylab = "Human liver clusters",
        margins = c(5,12))
#cexRow = 3, cexCol = 3,
#Rowv = NA)

# Aizarani comparisons
# Read in woodchuck data
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_TLH_rpca_mito50_dropletQC_scClustViz.RData")
# Read in ortholog table
geneNameTable <- read.table("~/Dropbox/Zoe/scf_version/make_gtf/orthofinder_sc2/homologene/collectedOrthofinderPairings.tsv",
                            sep = "\t",
                            header = TRUE)
# Choose clustering resolution
Idents(scSeurat) <- "integrated_snn_res.0.6"
# Get averages from variable genes
woodchuckClusterAverages <- AverageExpression(scSeurat,
                                              assays = "SCT",
                                              slot = "scale.data")
# Use all clusters
woodchuckClusterAverages$SCT <- na.omit(t(scale(t(woodchuckClusterAverages$SCT))))
# Grab gene names
uniqueHier <- row.names(woodchuckClusterAverages$SCT)
# Replace these gene names with the human name
uniqueHier <- as.data.frame(uniqueHier)
# Bind with geneNameTable
newNames <- dplyr::left_join(uniqueHier, geneNameTable, by = "uniqueHier")
# If mikado protein column is blank, replace with uniqueHier (non-coding genes)
newNames$mikado_final_sc2_stringent_noMito_protein <- ifelse(is.na(newNames$mikado_final_sc2_stringent_noMito_protein), newNames$uniqueHier, newNames$mikado_final_sc2_stringent_noMito_protein)
# If human or mouse has NA, replace with mikado name column
newNames$humanOneToOne <- ifelse(is.na(newNames$humanOneToOne), newNames$mikado_final_sc2_stringent_noMito_protein, newNames$humanOneToOne)
# Fix format
woodchuckClusterAverages <- woodchuckClusterAverages$SCT
# Add gene names as human
row.names(woodchuckClusterAverages) <- newNames$humanOneToOne
# Make sure formatted as list of dataframes
woodchuckClusterAverages <- as.data.frame(woodchuckClusterAverages)
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
aizaraniAverages <- as.data.frame(aizaraniAverages$SCT)
# Order by row name
woodchuckAveragesOrdered <- woodchuckClusterAverages[order(row.names(woodchuckClusterAverages)),]
aizaraniAveragesOrdered <- aizaraniAverages[order(row.names(aizaraniAverages)),]
# Now find intersecting genes
aizaraniMatches <- intersect(row.names(woodchuckAveragesOrdered),
                             row.names(aizaraniAveragesOrdered))
# Look at how many genes matched
length(aizaraniMatches)
# Make new matrices with only matching gene names
aizaraniAveragesCor <- aizaraniAveragesOrdered[aizaraniMatches,]
woodchuckAveragesCor <- woodchuckAveragesOrdered[aizaraniMatches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(aizaraniAveragesCor))), seq(ncol(as.matrix(woodchuckAveragesCor)))),
            mapply(function(i, j) cor.test(as.matrix(aizaraniAveragesCor)[, i], as.matrix(woodchuckAveragesCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(aizaraniAveragesCor)),
               dimnames=list(colnames(as.matrix(aizaraniAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(aizaraniAveragesCor)),
                dimnames=list(colnames(as.matrix(aizaraniAveragesCor)), colnames(as.matrix(woodchuckAveragesCor))))
heatmap(cors,
        main = "Pearson correlation of woodchuck single-cell vs Aizarani et al. liver data",
        xlab = "Woodchuck single-cell clusters",
        ylab = "Human liver clusters",
        margins = c(5,12))
#cexRow = 3, cexCol = 3,
#Rowv = NA)
pearVal <- cor(humanAveragesCor, woodchuckAveragesCor, method = "pearson")
heatmap(pearVal,
        main = "Pearson of woodchuck single-nuc vs human single-nuc",
        xlab = "Woodchuck TLH clusters",
        ylab = "Human liver clusters",
        margins = c(5,12))
#cexRow = 3, cexCol = 3,
#Rowv = NA)
# Do Spearman
spearVal <- cor(humanAveragesCor, woodchuckAveragesCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of woodchuck single-nuc vs human single-nuc",
        xlab = "Woodchuck TLH clusters",
        ylab = "Human liver clusters",
        margins = c(5,12))
#cexRow = 3, cexCol = 3,
#Rowv = NA)


# Compare single-nuc to single-cell
# Read in single-nuc woodchuck data
load("~/Dropbox/Zoe/scf_version/analysis/multiome/healthy_integration/L192_L212_multiome_onlyRNA_rpca_dropletQC_0.5nf_scClustViz.RData")
singleNuc <- scSeurat
# Read in single-cell woodchuck data
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_TLH_rpca_mito50_dropletQC_scClustViz.RData")
singleCell <- scSeurat
# Choose clustering resolution
Idents(singleNuc) <- "integrated_snn_res.0.6"
Idents(singleCell) <- "integrated_snn_res.0.6"
# Get averages from variable genes
nucClusterAverages <- AverageExpression(singleNuc,
                                        assays = "SCT",
                                        slot = "scale.data")
scClusterAverages <- AverageExpression(singleCell,
                                       assays = "SCT",
                                       slot = "scale.data")
# Use all clusters
nucClusterAverages <- as.data.frame(na.omit(t(scale(t(nucClusterAverages$SCT)))))
scClusterAverages <- as.data.frame(na.omit(t(scale(t(scClusterAverages$SCT)))))
# Order by row name
nucAvgOrd <- nucClusterAverages[order(row.names(nucClusterAverages)),]
scAvgOrd <- scClusterAverages[order(row.names(scClusterAverages)),]
# Now find intersecting genes
matches <- intersect(row.names(nucAvgOrd),
                     row.names(scAvgOrd))
# Look at how many genes matched
length(matches)
# Make new matrices with only matching gene names
nucAvgCor <- nucAvgOrd[matches,]
scAvgCor <- scAvgOrd[matches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(nucAvgCor))), seq(ncol(as.matrix(scAvgCor)))),
            mapply(function(i, j) cor.test(as.matrix(nucAvgCor)[, i], as.matrix(scAvgCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(nucAvgCor)),
               dimnames=list(colnames(as.matrix(nucAvgCor)), colnames(as.matrix(scAvgCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(nucAvgCor)),
                dimnames=list(colnames(as.matrix(nucAvgCor)), colnames(as.matrix(scAvgCor))))
heatmap(cors,
        main = "Pearson correlation of single-nuc VS single-cell",
        xlab = "Single-cell clusters",
        ylab = "Single-nuc clusters")

# Compare woodchuck PBMCs with woodchuck TLH
# Read in data
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_TLH_rpca_mito50_dropletQC_scClustViz.RData")
tlh <- scSeurat
# If single-nuc
Idents(tlh) <- "integrated_snn_res.0.6"
# If single-cell
Idents(tlh) <- "integrated_snn_res.0.6"
load("~/Dropbox/Zoe/scf_version/analysis/integration/L212_L202_PBMC_rpca_scClustViz.RData")
pbmc <- scSeurat
Idents(pbmc) <- "integrated_snn_res.0.6"
# Get averages from variable genes
tlhClusterAverages <- AverageExpression(tlh,
                                       assays = "SCT",
                                       slot = "scale.data")
pbmcClusterAverages <- AverageExpression(pbmc,
                                        assays = "SCT",
                                        slot = "scale.data")
# Use all clusters
tlhClusterAverages <- as.data.frame(na.omit(t(scale(t(tlhClusterAverages$SCT)))))
pbmcClusterAverages <- as.data.frame(na.omit(t(scale(t(pbmcClusterAverages$SCT)))))
# Order by row name
tlhAvgOrd <- tlhClusterAverages[order(row.names(tlhClusterAverages)),]
pbmcAvgOrd <- pbmcClusterAverages[order(row.names(pbmcClusterAverages)),]
# Now find intersecting genes
matches <- intersect(row.names(tlhAvgOrd),
                     row.names(pbmcAvgOrd))
# Look at how many genes matched
length(matches)
# Make new matrices with only matching gene names
tlhAvgCor <- tlhAvgOrd[matches,]
pbmcAvgCor <- pbmcAvgOrd[matches,]
# Do Pearson
tmp <- with(expand.grid(seq(ncol(as.matrix(tlhAvgCor))), seq(ncol(as.matrix(pbmcAvgCor)))),
            mapply(function(i, j) cor.test(as.matrix(tlhAvgCor)[, i], as.matrix(pbmcAvgCor)[, j]),
                   Var1, Var2))
cors <- matrix(unlist(tmp['estimate', ]), nrow=ncol(as.matrix(tlhAvgCor)),
               dimnames=list(colnames(as.matrix(tlhAvgCor)), colnames(as.matrix(pbmcAvgCor))))
pvals <- matrix(unlist(tmp['p.value', ]), nrow=ncol(as.matrix(tlhAvgCor)),
                dimnames=list(colnames(as.matrix(tlhAvgCor)), colnames(as.matrix(pbmcAvgCor))))
heatmap(cors,
        main = "Pearson correlation of woodchuck TLH VS PBMCs",
        xlab = "PBMC clusters",
        ylab = "scRNA-seq TLH clusters")
# Do spearman
spearVal <- cor(tlhAvgCor, pbmcAvgCor, method = "spearman")
heatmap(spearVal,
        main = "Spearman of woodchuck TLH vs woodchuck PBMCs",
        xlab = "Woodchuck PBMC clusters",
        ylab = "Woodchuck TLH clusters",
        margins = c(5,12))
