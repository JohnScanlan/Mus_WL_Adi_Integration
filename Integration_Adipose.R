library(Seurat)
library(SeuratDisk)
library(glmpca)
library(leiden)
library(leidenalg)
library(reticulate)
library(dplyr)

library(reticulate)
use_virtualenv("C:/Users/crtuser/leidenalg-env", required = TRUE)
leidenalg <- reticulate::import("leidenalg")

hasty <- readRDS("~/PhD/Project/Data/Hasty CiteSeq/GSE182233_IntegratedData.rds/GSE182233_IntegratedData.rds")
hasty <- subset(hasty, subset = orig.ident != 'WC')
cr_epi <- readRDS("~/PhD/Project/Data/Public_CR_scRNA/df.rds")
wl_epi <- LoadH5Seurat("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\MiceWeightLossCohort\\MOUSE RNA-seq\\preprocessed_epi.h5seurat")

hasty <- subset(hasty, idents = c('gd T Cells', 'ILC2s', 'Th2 CD4+ T Cells', 'NKT Cells', 'NK Cells', 'Cycling CD8+ T Cells', 'Central Memory CD8+ T Cells', 'Effector Memory CD8+ T Cells', 'Tregs', 'Th1 CD4+ T Cells'))
cr_epi <- subset(cr_epi, idents = c('21', '18', '16', '14', '11', '2', '5', '20'))
wl_epi <- subset(wl_epi, idents = c('yd_T', 'Activated_T', 'Gata3', 'T_regs', 'Tcf7_Cd4_T', 'Cd8_T', 'Cd27_NKs', 'Cd11b_NKs'))

hasty <- SCTransform(hasty, assay = 'RNA')
cr_epi <- SCTransform(cr_epi, assay = 'RNA')
wl_epi <- SCTransform(wl_epi, assay = 'RNA')

wl_epi@meta.data$source <- 'wl_epi'
cr_epi@meta.data$source <- 'cr_epi'
hasty@meta.data$source <- 'hasty'

VariableFeatures(hasty[["SCT"]]) <- rownames(hasty[["SCT"]]@scale.data)
VariableFeatures(cr_epi[["SCT"]]) <- rownames(cr_epi[["SCT"]]@scale.data)
VariableFeatures(wl_epi[["SCT"]]) <- rownames(wl_epi[["SCT"]]@scale.data)

# Create a list of Seurat objects
seurat_list <- list(hasty = hasty, cr_epi = cr_epi, wl_epi = wl_epi)

features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)

df <- IntegrateData(anchorset = immune.anchors)

DefaultAssay(df) <- "integrated"

DefaultAssay(object = df) <- "RNA"
df <- NormalizeData(object = df, assay = 'RNA')
df <- ScaleData(df, assay = 'RNA')

# Your existing code for preparing the GLMPCA data
sct_data <- df@assays$RNA$data
non_zero_rows <- rowSums(sct_data) > 0
sct_data_non_zero <- sct_data[non_zero_rows, ]

# Run GLMPCA on the subsetted data
glmpca_results <- glmpca(Y = sct_data_non_zero, L = 30, minibatch = 'stochastic')

colnames(glmpca_results$factors) <- paste0('glmpca', 1:30)

df[['glmpca']] <- CreateDimReducObject(embeddings = as.matrix(glmpca_results$factors), loadings = as.matrix(glmpca_results$loadings), assay = DefaultAssay(df), key = 'glmpca_')


#df <- ScaleData(df, verbose = FALSE)
#df <- RunGLMPCA(df, reduction.name = 'glmpca')
ElbowPlot(df, ndims = 1:30)
df <- RunUMAP(df, dims =  1:30, reduction = 'glmpca')
df <- FindNeighbors(df, reduction = 'glmpca')
df <- FindClusters(df, resolution = 0.7, algorithm = 4, method = 'igraph')

DimPlot(df, label = T, split.by = 'orig.ident') + NoLegend()

DotPlot(df, features = 'Tnfaip3', group.by = 'orig.ident')

#markers <- FindAllMarkers(df, logfc.threshold = 0.5)

FeaturePlot(df, features = 'Tnfaip3', label = T, pt.size = 1.3)

DefaultAssay(df) <- "RNA"
df <- NormalizeData(df, assay = 'RNA')
df <- ScaleData(df, assay = 'RNA')

# Define a mapping of original values to new values
condition_mapping <- list(
  "Lean" = c("epi_SFD"),
  "Obese" = c("ob", "epi_HFD"),
  "WL" = c("cr", "epi_WL")
)

# Add the 'Condition' column to your metdf
df@meta.data <- df@meta.data %>%
  mutate(Condition = case_when(
    orig.ident %in% condition_mapping[["Lean"]] ~ "Lean",
    orig.ident %in% condition_mapping[["Obese"]] ~ "Obese",
    orig.ident %in% condition_mapping[["WL"]] ~ "WL",
    TRUE ~ as.character(orig.ident)
  ))


#saveRDS(df, 'C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Integrated_Mice_Adi_WL\\louvain_integration.rds')


# Load necessary libraries
library(Seurat)
library(glmpca)

# Assume you have Seurat objects: hasty, cr_epi, wl_epi
# Add a source label to each
hasty$source <- 'hasty'
cr_epi$source <- 'cr_epi'
wl_epi$source <- 'wl_epi'

# Create a list of Seurat objects
seurat_list <- list(hasty = hasty, cr_epi = cr_epi, wl_epi = wl_epi)

# Integration Steps
# Normalize and identify variable features for each dataset
for (name in names(seurat_list)) {
  seurat_list[[name]] <- SCTransform(seurat_list[[name]], verbose = TRUE, vars.to.regress = 'percent.mt', method = 'glmGamPoi')
  seurat_list[[name]] <- RunPCA(seurat_list[[name]], features = VariableFeatures(seurat_list[[name]]), assay = "SCT", npcs = 50, verbose = TRUE)
}

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)

# Preparing objects for integration
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

# Finding integration anchors
immune.anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features, reduction = "rpca")

# Integrating the data
cells.combined <- IntegrateData(anchorset = immune.anchors, new.assay.name = "integrated", normalization.method = 'SCT', verbose = TRUE)
DefaultAssay(cells.combined) <- "integrated"

# Prepare Data for GLMPCA
# Extract integrated data
integrated_data <- cells.combined@assays$integrated@data

# Filter out rows with all zeros
non_zero_rows <- rowSums(integrated_data) > 0
integrated_data_non_zero <- integrated_data[non_zero_rows, ]

# Run GLMPCA
glmpca_results <- glmpca(Y = integrated_data_non_zero, L = 30, minibatch = 'stochastic')
colnames(glmpca_results$factors) <- paste0('glmpca', 1:30)

# Add GLMPCA Results to Seurat Object
cells.combined[['glmpca']] <- CreateDimReducObject(
  embeddings = as.matrix(glmpca_results$factors),
  loadings = as.matrix(glmpca_results$loadings),
  assay = DefaultAssay(cells.combined),
  key = 'glmpca_'
)

# cells.combined now contains the integrated data and GLMPCA results
