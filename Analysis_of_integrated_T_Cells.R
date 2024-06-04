library(BITFAM)
library(Seurat)
library(scCustomize)
library(ggvolc)
library(data.table)
library(dplyr)
library(reshape2)  # or library(data.table) if you're using that for dcast
library(ggplot2)
library(scales)
library(tidyverse)
library(scales)
library(RColorBrewer)
library(clusterProfiler)
library(data.table)
library(org.Mm.eg.db)
library(ggplot2)
library(ggnewscale)
library(GOSemSim)
library(DOSE)
library(enrichplot)

#Load in the BITFAM (SC) and RNA-seq file
#sc <- readRDS("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Integrated_Mice_Adi_WL\\bitfam_integrated_adipose_T.rds")
df <- readRDS("C:\\Users\\crtuser\\OneDrive - TCDUD.onmicrosoft.com\\Documents\\PhD\\Project\\Data\\Integrated_Mice_Adi_WL\\louvain_integration.rds")
df@meta.data$cell_type <- Idents(df)
Idents(df) <- df@meta.data$cell_type
# Assuming 'df' is your Seurat object
nks <- subset(df, subset = cell_type == 'Cd11b NK Cells')
# If df is a data frame, ensure the syntax is correct
# Extract metadata


FeaturePlot(df, features = 'Kcnip3', pt.size = 1.2, label = T, order = F)
VlnPlot(df, features = 'Cpt1a', pt.size = 0, group.by = 'Condition') + NoLegend()
DotPlot(df, features = c('Stat1', 'Ifng', 'Ifngr1', 'Slc3a2', 'Tnfrsf1a', 'Tnfrsf1b', 'Kcnip3')) + RotatedAxis()
ggsave('C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Writing/Figures/Tnfaip3_Vln.png', f, width = 8, height = 4, dpi = 400)

View(df@meta.data)

a20 <- FindMarkers(df, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Tnfaip3_exp', max.cells.per.ident = 1000)
View(a20)

sc[['umap']] <- df[['umap']]

# Pathway Analysis
# TNFAIP3 Positive Cells

#Significant Clusters
# Assuming 'df' is a Seurat object already created and preprocessed

# Define the vector of cell types
cells <- c("Cd27 NKs Cells", "ZFP36+ ILCs", "Cd11b NK Cells", "Gamma Delta T Cells", "Cd8 T Cells", "Activated CD4 T Cells",
           "Naive T Cells", "PD1+ CD4 T Cells", "SRM+ NK Cells", "B Cells", "TXNIP+ ILCs", "Proliferating T Cells", "T-Regs", "Myeloid", "Granulocytes")

ms <- FindAllMarkers(df, logfc.threshold = 0.75)
View(ms)

# Find Markers - NK Cells
# Alternative subsetting using WhichCells and subset


m <- FindMarkers(subset(df, idents = 'Cd11b NK Cells'), group.by = 'Condition', ident.1 = 'WL', ident.2 = 'Obese')
nks <- subset(df, idents = 'Cd11b NK Cells')
# Create an empty data frame to store the combined results
combined_results <- data.frame()

# Iterate over each cell type
for (cell_type in cells) {
  # Check if the cell type exists in the Seurat object
  if(cell_type %in% levels(Idents(df))) {
    # Perform the analysis for the current cell type
    m.markers <- FindMarkers(subset(df, idents = cell_type), group.by = 'Condition', ident.1 = 'WL', ident.2 = 'Obese', logfc.threshold = 0.0)
    
    # Add the cell type as a new column
    m.markers$Cell_Type <- cell_type
    
    # Append the current m.markers data frame to the combined results
    combined_results <- rbind(combined_results, m.markers)
  } else {
    warning(paste("Cell type", cell_type, "not found in Seurat object. Skipping..."))
  }
}

# Print the combined results
# Replace View with write.csv if you are in a non-interactive environment
View(combined_results)

# Assuming combined_results is your original data frame with all the data.
# Extract only rows where the row names contain 'Tnfaip3'
tnfaip3_rows <- grep("Tnfaip3", rownames(combined_results), value = TRUE)

# Use the extracted row names to subset the original data frame
tnfaip3_df <- combined_results[tnfaip3_rows, ]

# Remove the unwanted cell types: granulocytes, myeloid cells, and B cells
tnfaip3_df <- tnfaip3_df[!tnfaip3_df$Cell_Type %in% c("Granulocytes", "Myeloid", "B Cells"),]

# Now we can plot using ggplot2
library(ggplot2)

# Define custom labels for the significance levels
tnfaip3_df$stars <- ifelse(tnfaip3_df$p_val_adj <= 0.001, '***', 
                           ifelse(tnfaip3_df$p_val_adj <= 0.01, '**', 
                                  ifelse(tnfaip3_df$p_val_adj <= 0.05, '*', '')))

# Define custom labels for the fill based on the sign of avg_log2FC
tnfaip3_df$Label <- ifelse(tnfaip3_df$avg_log2FC > 0, 'Increased in WL', 'Increased in Obese')

# Calculate y positions for the significance stars, make sure it is above the highest bar
max_value <- max(tnfaip3_df$avg_log2FC, na.rm = TRUE)
star_height <- max_value + (max_value * 0.1)  # Increase by 10% of max value

# Add a new column for the y position of the stars
tnfaip3_df$star_y <- ifelse(!is.na(tnfaip3_df$stars), star_height, NA)

p <- ggplot(tnfaip3_df, aes(x = reorder(Cell_Type, avg_log2FC), y = avg_log2FC, fill = Label)) +
  geom_col(color = "black", linewidth = 0.25) + # Add border to bars
  scale_fill_manual(values = c("Increased in WL" = "#56B4E9", "Increased in Obese" = "red")) +
  labs(title = "Tnfaip3 Expression - WL vs Obese",
       y = "Average Log2 Fold Change") +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black", size = 12, face = "bold"),
    axis.text.y = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  geom_text(aes(label = stars, y = star_y), size = 5, color = "black", vjust = 0)

# Show the plot
print(p)


#### BULK
df <- subset(df, idents = c('Myeloid', 'B Cells', 'Granulocytes'), invert = T)

# Filter the 'bulk' dataframe for the specified genes
selected_genes <- bulk[rownames(bulk) %in% c("Fabp4", "Tigit", "Lars2", "Nr4a2", "Sqstm1", "Nme2", "Ube2v1", 'Fos', ),]

# Add a column for gene names if it does not already exist
selected_genes$Gene <- rownames(selected_genes)

# Create a column for Regulation based on avg_log2FC for coloring
selected_genes$Regulation <- ifelse(selected_genes$avg_log2FC > 0, 'Upregulated', 'Downregulated')

# Now we plot using ggplot2
library(ggplot2)

# Create the bar chart with enhanced aesthetics
p <- ggplot(selected_genes, aes(x = Gene, y = avg_log2FC, fill = Regulation)) +
  geom_bar(stat = "identity") + # This will create the bars
  scale_fill_manual(values = c('Upregulated' = '#00BA38', 'Downregulated' = '#F8766D')) + # Set colors for up/down
  labs(title = "Log2 Fold Change of Selected Genes",
       x = "Gene",
       y = "Average Log2 Fold Change") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add a dashed line at y = 0
  coord_flip() # Flip the coordinates for horizontal bars

# Show the plot
print(p)

#Tnfaip3 - Pos vs Neg
# Non Canonical
poscells <- WhichCells(df, expression = Tnfaip3 > 0.75)
df@meta.data$Tnfaip3_exp <- ifelse(colnames(df) %in% poscells, "Pos", "Neg")

gene_markers <- FindMarkers(df, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Tnfaip3_exp', logfc.threshold = 0)

genes_of_interest <- c("Nfkb1",    # NF-??B1 p50 subunit, part of the NF-??B complex
                       "Rela",     
                       "Ikbkb",    
                       "Ikbkg",    
                       "Tnf",  'Hif1a',  
                       "Il1b",    
                       "Ripk1",  
                       "Traf6",   
                       "Rel", 'Birc3', 'Map3k8', 'Ikbke',
                       "Tnfrsf1a",
                       "Tnfrsf1b", 
                       "Nfkbia", "Nfkb1", 'Relb', 'Nfkb2', 
                       'Chuk', 'Traf6', 'Traf5', 'Traf3', 
                       'Traf4', 'Il17ra', 'Il17rc', 'Tradd',
                       'Tax1bp1', 'Tab1', 'Tab2', 'Tab3', 'Tnip1', 'Tnip2', 'Tnip3', 'Il6', 'Nfkbiz', 'Nfkbid', 'Ikbke')  

gene_markers$X <- rownames(gene_markers)
#attention <- subset(gene_markers, X %in% genes_of_interest)
attention <- subset(gene_markers, X == genes_of_interest)

gene_markers$pvalue <- gene_markers$p_val_adj
gene_markers$log2FoldChange <- gene_markers$avg_log2FC
gene_markers$genes <- gene_markers$X

genes_to_remove <- c('Tnfaip3')

# Filter the data frame
gene_markers <- gene_markers[!gene_markers$X %in% genes_to_remove, ]
attention_genes <- gene_markers[gene_markers$X %in% genes_of_interest, ]

#attention <- subset(gene_markers, X %in% genes_of_interest)
#attention <- subset(gene_markers, X == genes_of_interest)

ggvolc(gene_markers, attention_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="NF\u03BAB Genes - Tnfaip3 Positive vs Negative T Cells")

# Canonical
genes_of_interest <- c("Nfkb1",    # NF-??B1 p50 subunit, part of the NF-??B complex
                                 "Rela",     
                                 "Ikbkb",    
                                 "Ikbkg",    
                                 "Tnf",   
                                  "Il1b",    
                                 "Ripk1",  
                                 "Traf6", 'Ldha',   
                                 "Rel",
                                "Tnfrsf1a",
                                "Tnfrsf1b", 
                                "Nfkbia", "Nfkb1", 'Relb', 'Nfkb2', 
                                'Chuk', 'Traf6', 'Traf5', 'Traf3', 
                                'Traf4', 'Il17ra', 'Il17rc', 'Tradd',
                                'Tax1bp1', 'Tab1', 'Tab2', 'Tab3', 'Tnip1', 'Tnip2', 'Tnip3', 'Il6', 'Nfkbiz', 'Nfkbid', 'Ikbke')      

gene_markers$X <- rownames(gene_markers)
#attention <- subset(gene_markers, X %in% genes_of_interest)
attention <- subset(gene_markers, X == genes_of_interest)

gene_markers$pvalue <- gene_markers$p_val_adj
gene_markers$log2FoldChange <- gene_markers$avg_log2FC
gene_markers$genes <- gene_markers$X

genes_to_remove <- c('Tnfaip3')

# Filter the data frame
gene_markers <- gene_markers[!gene_markers$X %in% genes_to_remove, ]
attention_genes <- gene_markers[gene_markers$X %in% genes_of_interest, ]

#attention <- subset(gene_markers, X %in% genes_of_interest)
#attention <- subset(gene_markers, X == genes_of_interest)

ggvolc(gene_markers, attention_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="NF\u03BAB Genes - Tnfaip3 Positive vs Negative T Cells")

#Top Genes
# Sort the gene_markers by adjusted p-value
#Negatively and positively  associated genes with Tnfaip3
# Sort gene_markers by positive log fold change and select the top 20
positive_genes <- gene_markers[gene_markers$avg_log2FC > 0, ]
top_positive_genes <- head(positive_genes[order(-positive_genes$avg_log2FC), ], 80)

# Sort gene_markers by negative log fold change and select the top 10
negative_genes <- gene_markers[gene_markers$avg_log2FC < 0, ]
top_negative_genes <- head(negative_genes[order(negative_genes$avg_log2FC), ], 15)

# Combine the top positive and top negative genes
combined_genes <- rbind(top_positive_genes, top_negative_genes)

# Remove any genes you don't want to include
genes_to_remove <- c('Tnfaip3')
combined_genes <- combined_genes[!rownames(combined_genes) %in% genes_to_remove, ]

# Update the combined_genes dataframe with the required columns for plotting
combined_genes$pvalue <- combined_genes$p_val_adj
combined_genes$log2FoldChange <- combined_genes$avg_log2FC
combined_genes$genes <- rownames(combined_genes)

# Now you can create your plot using the combined_genes dataframe
# Your existing ggvolc plotting code would be modified to use combined_genes
ggvolc(gene_markers, combined_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="Top Genes - Tnfaip3 Positive and Negative Cells")

#Just WL:
# Assuming df is your SingleCellExperiment object and has a metadata column named 'Condition'
wl_cells <- df[, df@meta.data$Condition == 'WL']

# Now, run FindMarkers on the subsetted data
gene_markers <- FindMarkers(wl_cells, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Tnfaip3_exp', logfc.threshold = 0)

# Assuming gene_markers is the output from FindMarkers
gene_markers$genes <- rownames(gene_markers)  # Add a new column with row names

# Separate into positive and negative associated genes
positive_genes <- gene_markers[gene_markers$avg_log2FC > 0, ]
negative_genes <- gene_markers[gene_markers$avg_log2FC < 0, ]

# Sort and select top 10 positive and negative genes
top_positive_genes <- head(positive_genes[order(-positive_genes$avg_log2FC), ], 10)
top_negative_genes <- head(negative_genes[order(negative_genes$avg_log2FC), ], 10)

# Combine the top genes from both positive and negative sets
top_genes <- rbind(top_positive_genes, top_negative_genes)

# Remove 'Tnfaip3' from the top genes if it's present
genes_to_remove <- c('Tnfaip3')
top_genes <- top_genes[!top_genes$genes %in% genes_to_remove, ]

# Adjust the column names to match ggvolcano's expected input
top_genes$pvalue <- top_genes$p_val_adj
top_genes$log2FoldChange <- top_genes$avg_log2FC
# Plot the results using ggvolcano (assuming ggvolc is a typo for ggvolcano)
# Note: ggvolcano might not be a standard function, ensure you have the correct package or function for plotting
ggvolc(gene_markers, top_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="WL - Tnfaip3 Positive vs Negative T Cells")

## Obesity
# Correctly subset the Seurat object for 'Obese' condition
obese <- df[, df@meta.data$Condition == 'Obese']

# Now, run FindMarkers on the subsetted data
gene_markers <- FindMarkers(obese, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Tnfaip3_exp', logfc.threshold = 0)
gene_markers$pvalue <- gene_markers$p_val_adj
gene_markers$log2FoldChange <- gene_markers$avg_log2FC
# Assuming gene_markers is the output from FindMarkers
gene_markers$genes <- rownames(gene_markers)  # Add a new column with row names

# Separate into positive and negative associated genes
positive_genes <- gene_markers[gene_markers$avg_log2FC > 0, ]
negative_genes <- gene_markers[gene_markers$avg_log2FC < 0, ]

# Sort and select top 10 positive and negative genes
top_positive_genes <- head(positive_genes[order(-positive_genes$avg_log2FC), ], 20)
top_negative_genes <- head(negative_genes[order(negative_genes$avg_log2FC), ], 10)

# Combine the top genes from both positive and negative sets
top_genes <- rbind(top_positive_genes, top_negative_genes)

# Remove 'Tnfaip3' from the top genes if it's present
genes_to_remove <- c('Tnfaip3')
top_genes <- top_genes[!top_genes$genes %in% genes_to_remove, ]
gene_markers <- gene_markers[!gene_markers$genes %in% genes_to_remove, ]

# Adjust the column names to match ggvolcano's expected input
top_genes$pvalue <- top_genes$p_val_adj
top_genes$log2FoldChange <- top_genes$avg_log2FC
# Plot the results using ggvolcano (assuming ggvolc is a typo for ggvolcano)
# Note: ggvolcano might not be a standard function, ensure you have the correct package or function for plotting
ggvolc(gene_markers, top_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="Obese - Tnfaip3 Positive vs Negative T Cells")


#Cell Proportions
md <- df@meta.data %>% as.data.table
# the resulting md object has one "row" per cell

## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
md[, .N, by = c("Condition", "Cell")]
View(md)

md <- md[, .N, by = c("Condition", "Cell")] %>% dcast(., Condition ~ Cell, value.var = "N")
View(md)

ggplot(md, aes(x = Cell, y = Condition)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Pastel1")

# Assuming you've already processed the data using dcast
md <- md[, .N, by = c("Condition", "Cell")] %>% dcast(., Condition ~ Cell, value.var = "N")

# Assuming md is your data frame after some transformations
md <- md %>% 
  mutate_at(vars(-Condition), ~ . / sum(.)) %>%
  gather(Cell, proportion, -Condition)

# Determine the number of unique Cells
number_of_cells <- length(unique(md$Cell))

# Create a combined palette
# If 'number_of_cells' is large, you may want to generate colors differently
# For demonstration, combining 'Set1', 'Set2', and 'Set3'
all_palettes <- c(brewer.pal(9, "Set1"), brewer.pal(9, "Set2"), brewer.pal(12, "Set3"))
combined_palette <- all_palettes[1:min(number_of_cells, length(all_palettes))]

# Create the stacked bar plot using ggplot2 with the combined palette
plots <- ggplot(md, aes(x = Condition, y = proportion, fill = Cell)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_y_continuous(labels = percent_format(scale = 100)) +
  scale_fill_manual(values = combined_palette) + # Using manual scale for colors
  labs(title = "Proportion of Cells in Each Experimental Group",
       y = "Proportion (%)",
       x = "Experimental Group",
       fill = "Cell Group") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "right", # or "bottom" if you have many groups
    legend.key.size = unit(0.5, "cm"), # Adjust legend key size
    plot.margin = margin(10, 10, 10, 10) # Adjust plot margins (top, right, bottom, left)
  )

# Save the plot to a file (good for publication or sharing)
ggsave("my_plot.png", plots, width = 11, height = 8, dpi = 300)

# Print the plot to R's current graphics device
print(plots)


### BITFAM
View(sc@meta.data)

gene_markers <-FindMarkers(sc, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Tnfaip3_exp', logfc.threshold = 0)

# Assuming gene_markers is the output from FindMarkers
gene_markers$genes <- rownames(gene_markers)  # Add a new column with row names

# Separate into positive and negative associated genes
positive_genes <- gene_markers[gene_markers$avg_log2FC > 0, ]
negative_genes <- gene_markers[gene_markers$avg_log2FC < 0, ]

# Sort and select top 10 positive and negative genes
top_positive_genes <- head(positive_genes[order(-positive_genes$avg_log2FC), ], 10)
top_negative_genes <- head(negative_genes[order(negative_genes$avg_log2FC), ], 10)

# Combine the top genes from both positive and negative sets
top_genes <- rbind(top_positive_genes, top_negative_genes)

# Remove 'Tnfaip3' from the top genes if it's present
genes_to_remove <- c('Tnfaip3')
top_genes <- top_genes[!top_genes$genes %in% genes_to_remove, ]

# Adjust the column names to match ggvolcano's expected input
top_genes$pvalue <- top_genes$p_val_adj
top_genes$log2FoldChange <- top_genes$avg_log2FC
gene_markers$pvalue <- gene_markers$p_val_adj
gene_markers$log2FoldChange <- gene_markers$avg_log2FC
# Plot the results using ggvolcano (assuming ggvolc is a typo for ggvolcano)
# Note: ggvolcano might not be a standard function, ensure you have the correct package or function for plotting
ggvolc(gene_markers, top_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="BITFAM - Tnfaip3 Positive vs Negative T Cells")

#Pathway Analysis
performPathwayAnalysis <- function(df, ident1, ident2, groupby, orgDb) {
  # Find markers between two identities
  markers <- FindMarkers(df, 
                         ident.1 = ident1, 
                         ident.2 = ident2,
                         group.by = groupby,
                         logfc.threshold = 0.1)
  
  # View the markers
  View(markers)
  
  # Add gene names to the markers data frame
  markers$genes <- rownames(markers)

  # Subset for significant markers
  pos_markers <- subset(markers, p_val_adj < 0.1 & avg_log2FC > 0.1)
  neg_markers <- subset(markers, p_val_adj < 0.1 & avg_log2FC < -0.1)

  # Extract gene names for pathway analysis
  gene_names <- pos_markers$genes
  print(gene_names)
  # Perform pathway enrichment analysis
  enrich_result <- enrichGO(gene = gene_names,
                            keyType = 'SYMBOL',
                            OrgDb = orgDb, pvalueCutoff = 0.1)
  print(enrich_result)
  # Visualize enriched pathways with a dot plot
  dot_plot <- dotplot(enrich_result, showCategory = 10) +
    ggtitle(paste(ident1, 'vs', ident2, 'Enriched Pathways'))
  
  # Return the dot plot and the enrichment result
  list(dot_plot = dot_plot, enrich_result = enrich_result)
}

# Example usage:
# Assuming 'df' is your data frame, 'gene_info' contains gene annotations,
# and 'org.Mm.eg.db' is the gene annotations package for mice.
# Replace 'ly6c_pos', 'ly6c_neg', and 'ly6c_status' with the actual values you want to use.

result <- performPathwayAnalysis(df = df, 
                                 ident1 = 'Obese', 
                                 ident2 = 'WL', 
                                 groupby = 'Condition', 
                                 orgDb = org.Mm.eg.db)

# Differential Expression - Ly6c+ vs Ly6c-
gene_markers <- FindMarkers(df, 
                            ident.1 = 'Pos', 
                            ident.2 = 'Neg',
                            group.by = 'Tnfaip3_exp',
                            logfc.threshold = 0.25)
View(gene_markers)

gene_names_mouse <- rownames(gene_markers)

#DataFrame is named 'df', use this to access files
gene_markers_pos <- subset(gene_markers, p_val_adj < 0.05 & avg_log2FC > 0.1)
gene_markers_neg <- subset(gene_markers, p_val_adj < 0.05 & avg_log2FC < -0.1)

log2_fold_changes_mouse <- gene_markers$avg_log2FC
p_values_mouse <- gene_markers$p_val_adj

# Load gene annotations package for mice
gene_annot <- org.Mm.eg.db

# Perform pathway enrichment analysis
enrich_result_mouse <- enrichGO(gene = rownames(gene_markers_neg),
                                keyType = 'SYMBOL',
                                OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05)

# Visualize enriched pathways with a dot plot
dotplot(enrich_result_mouse, showCategory = 10) +
  ggtitle(expression(paste('Tnfaip3 Negative Gene Associations - Enriched Pathways')))

heatplot(enrich_result_mouse, showCategory=10, foldChange = gene_markers_neg$avg_log2FC)




###### CD4 T CELLS
gene_markers <- FindMarkers(cd4, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Tnfaip3_exp', logfc.threshold = 0)
gene_markers$pvalue <- gene_markers$p_val_adj
gene_markers$log2FoldChange <- gene_markers$avg_log2FC
# Assuming gene_markers is the output from FindMarkers
gene_markers$genes <- rownames(gene_markers)  # Add a new column with row names

# Separate into positive and negative associated genes
positive_genes <- gene_markers[gene_markers$avg_log2FC > 0, ]
negative_genes <- gene_markers[gene_markers$avg_log2FC < 0, ]

# Sort and select top 10 positive and negative genes
top_positive_genes <- head(positive_genes[order(-positive_genes$avg_log2FC), ], 10)
top_negative_genes <- head(negative_genes[order(negative_genes$avg_log2FC), ], 10)

# Combine the top genes from both positive and negative sets
top_genes <- rbind(top_positive_genes, top_negative_genes)

# Remove 'Tnfaip3' from the top genes if it's present
genes_to_remove <- c('Tnfaip3')
top_genes <- top_genes[!top_genes$genes %in% genes_to_remove, ]
gene_markers <- gene_markers[!gene_markers$genes %in% genes_to_remove, ]

# Adjust the column names to match ggvolcano's expected input
top_genes$pvalue <- top_genes$p_val_adj
top_genes$log2FoldChange <- top_genes$avg_log2FC
# Plot the results using ggvolcano (assuming ggvolc is a typo for ggvolcano)
# Note: ggvolcano might not be a standard function, ensure you have the correct package or function for plotting
ggvolc(gene_markers, top_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="CD4 T Cells - Tnfaip3 Positive vs Negative T Cells")


data <- FetchData(df, vars = c('orig.ident', 'Tnfaip3', 'Nfkb1', 'Rela', 'Rel', 'Nfkb2', 'Relb', 'Condition', 'clusters'))

pseudobulk_data <- data %>%
  group_by(Condition, clusters) %>%
  summarise(across(Tnfaip3:Nfkb2, mean, na.rm = TRUE)) # Replace `mean` with `sum` if preferred

# View the pseudobulk data
print(pseudobulk_data)

library(dplyr)
library(ggplot2)

library(dplyr)
library(tidyr)
library(ggplot2)

# Assuming 'data' is your initial dataframe
# Transform data to long format
pseudobulk_long <- data %>%
  pivot_longer(cols = Tnfaip3:Nfkb2, names_to = "Gene", values_to = "Expression") %>%
  mutate(
    Condition = factor(Condition, levels = c("Lean", "Obese", "WL")),
    clusters = factor(clusters) # Convert clusters to factor if not already
  )

# This is a simple approach for demonstration. Adjust as needed for your specific ordering requirements.
pseudobulk_long <- pseudobulk_long %>%
  arrange(Condition, clusters) %>%
  mutate(clusters = factor(clusters, levels = unique(clusters)))

ggplot(pseudobulk_long, aes(x = clusters, y = Expression, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Condition, ncol = 1, scales = "free_x") +
  theme_minimal() +
  labs(title = "Gene Expression by Cell Type and Condition",
       x = "Cell Type",
       y = "Expression") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  )


# Split data by 'Condition'
data_by_condition <- split(data, data$Condition)

# Names of genes for selecting columns
gene_names <- c("Tnfaip3", "Nfkb1", "Rela", "Rel", "Nfkb2", "Relb")

for (condition_name in names(data_by_condition)) {
  # Select gene columns and calculate correlation matrix
  gene_data <- dplyr::select(data_by_condition[[condition_name]], gene_names)
  corr_matrix <- cor(gene_data, use = "pairwise.complete.obs", method = "pearson")
  
  # Remove self-correlations by setting diagonal to NA
  diag(corr_matrix) <- NA
  
  # Plot the correlation matrix
  corrplot::corrplot(corr_matrix, method = "color", type = "upper", order = "hclust",
                     tl.col = "black", tl.srt = 45, addCoef.col = "black",
                     title = paste("Correlation -", condition_name))
}

Stacked_VlnPlot(seurat_object = df, features = c("Tnfaip3", "Nfkb1"),
                x_lab_rotate = TRUE,
                split.by = 'Condition',
                plot_legend = T) 

library(ggplot2)

# Generate the plot
plot_object <- Stacked_VlnPlot(seurat_object = df, features = c("Tnfaip3", "Nfkb1"),
                               x_lab_rotate = TRUE,
                               split.by = 'Condition',
                               plot_legend = TRUE, 
                               colors_use = c('black', 'red', 'steelblue'))

# Inspect the plot structure
plot_object
ggsave("C:/Users/crtuser/OneDrive - TCDUD.onmicrosoft.com/Documents/PhD/Project/Writing/Figures/Tnfaip3_Nfkb1_Vln.png", plot_object, width = 8, height = 4, dpi = 400)


# Adjust the color/fill scale
plot_object + scale_fill_manual(values = c("Condition1" = "black", "Condition2" = "red", "Condition3" = "steelblue"))
# Or if it uses color instead of fill
plot_object + scale_color_manual(values = c("Condition1" = "black", "Condition2" = "red", "Condition3" = "steelblue"))

genes <- c('ZFP36','CCL1', "CCL11", "CCL15", "CCL17", "CCL19", "CCL2", "CCL20", "CCL21", "CCL22", "CCL23", 
           "CCL24", "CCL3", "CCL4", "CCL5", "CCL6", "CCR1", "CCRL2", "CX3CL1", "CX3CR1", "CXCL1", 
           "CXCL10", "CXCL11", "CXCL12", "CXCL16", "CXCL2", "CCL28", "CXCL6", "CXCL8", "CXCL9", 
           "CXCR4", "CXCR4", "CXCR5", "EBI3", "ICOS", "IFNAB", "IFNAR1", "IFNB1", "IFNE", "IFNG", 
           "IFNGR1", "IFNGR2", "IFNK", "IFNL1", "IIGP1", "IL10", "IL10RA", "IL10RB", "IL11", "IL12A", 
           "IL12B", "IL12RB1", "IL13", "IL13RA1", "IL15", "IL15RA", "IL1A", "IL1B", "IL1R1", "IL1R2", 
           "IL1RAP", "IL1RN", "IL2", "IL21", "IL23A", "IL27", "IL32", "IL4", "IL4I1", "IL4R", "IL6", 
           "IL6ST", "IL7R", "IL8", "IL9", "LIFR", "LTA", "LTB", "TFF3", "TNF", "TNFRSF10B", "TNFSF10", 
           "TNFSF11", "TNFSF13B", "TNFSF14", "TNFSF4", "XCL1")

# Convert gene names to have only the first letter capitalized
formatted_genes <- sapply(genes, function(x) {
  paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2)), sep = "")
})

# Removing duplicates, as CXCR4 appeared twice
formatted_genes <- unique(formatted_genes)


# Correctly format gene names to have only the first letter capitalized
formatted_genes <- sapply(genes, function(x) {
  paste(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))), sep = "")
})

# Removing duplicates, as CXCR4 appeared twice
formatted_genes <- unique(formatted_genes)


formatted_genes

data <- FetchData(df, vars = c('orig.ident', 'Tnfaip3', 'Nfkb1', 'Rela', 'Rel', 'Nfkb2', 'Relb',formatted_genes,'Condition', 'clusters'))

# Assuming your data frame is named 'data'
# Filter out non-numerical columns
numerical_data <- data[, sapply(data, is.numeric)]

# Ensure 'Tnfaip3' is included if it was excluded by the numeric filter
# This step assumes 'Tnfaip3' is numeric; if not, this approach needs adjusting
if (!"Tnfaip3" %in% names(numerical_data)) {
  warning("'Tnfaip3' is not a numeric column. Please check your data.")
} else {
  # Calculate correlation of all numeric genes with Tnfaip3
  correlations <- sapply(numerical_data, function(x) cor(x, numerical_data$Tnfaip3, use = "complete.obs"))
  
  # Remove the correlation of Tnfaip3 with itself to avoid a perfect correlation entry
  correlations <- correlations[-which(names(correlations) == "Tnfaip3")]
  
  # Create a data frame to better visualize the gene names alongside their correlation with Tnfaip3
  correlation_df <- data.frame(Gene = names(correlations), CorrelationWithTnfaip3 = correlations)
  
  # Viewing the correlations sorted by their absolute values to see stronger correlations first
  correlation_df_sorted <- correlation_df[order(-abs(correlation_df$CorrelationWithTnfaip3)), ]
  
  print(correlation_df_sorted)
}

# Assuming df is your SingleCellExperiment or Seurat object
# And assuming 'df[['RNA']]@data' is a valid way to access the expression matrix

# Step 1: Extract the gene expression matrix
expression_matrix <- as.matrix(df[['RNA']]@data)

# Step 2: Extract Tnfaip3 expression vector
tnfaip3_expression <- expression_matrix["Tnfaip3", ]

# Step 3: Calculate correlations between Tnfaip3 and all genes
correlations <- apply(expression_matrix, 1, function(gene_expression) {
  cor(tnfaip3_expression, gene_expression, use = "pairwise.complete.obs")
})

# Step 4: Create a dataframe for better visualization and sorting of correlations
correlation_df <- data.frame(Gene = rownames(expression_matrix), CorrelationWithTnfaip3 = correlations)

# Optional: Sort by absolute correlation value to see the strongest correlations
correlation_df <- correlation_df[order(-abs(correlation_df$CorrelationWithTnfaip3)), ]

# View the sorted correlations
head(correlation_df)

###
# Placeholder for subsets - replace with your actual subsetting code
conditions <- c("Lean", "Obese", "WL") # Example condition names
condition_vectors <- list(Lean = lean_vector, Obese = obese_vector, WL = wl_vector) # Replace with actual vectors identifying columns for each condition

# Initialize a list to store correlation results for each condition
correlation_results <- list()

for (condition in conditions) {
  # Extract the subset of expression data for the current condition
  subset_expression_matrix <- expression_matrix[, condition_vectors[[condition]]]
  
  # Extract Tnfaip3 expression for the subset
  tnfaip3_expression_subset <- subset_expression_matrix["Tnfaip3", ]
  
  # Calculate correlations for the subset
  correlations_subset <- apply(subset_expression_matrix, 1, function(gene_expression) {
    cor(tnfaip3_expression_subset, gene_expression, use = "pairwise.complete.obs")
  })
  
  # Store correlations in the results list
  correlation_results[[condition]] <- correlations_subset
}

# Note: Ensure that 'lean_vector', 'obese_vector', and 'wl_vector' are correctly defined to index samples by condition
# Assuming the row names are the same across all conditions, use the row names from the first subset as gene names
gene_names <- rownames(subset_expression_matrix)

# Combine correlation results into a matrix, where each column represents a condition
correlation_matrix <- do.call(cbind, correlation_results)
rownames(correlation_matrix) <- gene_names
colnames(correlation_matrix) <- conditions

library(ggplot2)
library(ComplexHeatmap)

# Create the heatmap
Heatmap(correlation_matrix,
        name = "Correlation",
        top_annotation = HeatmapAnnotation(condition = colnames(correlation_matrix), col = list(condition = c(Lean = "green", Obese = "red", WL = "blue"))),
        left_annotation = rowAnnotation(Genes = anno_text(gene_names, rot = 45, offset = -1)),
        cluster_rows = TRUE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        column_title = "Correlation with TNFAIP3 by Condition",
        row_title = "Genes")



####-
library(Seurat)

# Subset Seurat object by condition
Lean <- subset(df, subset = Condition == 'Lean')
Obese <- subset(df, subset = Condition == 'Obese')
WL <- subset(df, subset = Condition == 'WL')

calculate_correlations <- function(seurat_object) {
  expression_matrix <- as.matrix(seurat_object[['RNA']]@data)
  tnfaip3_expression <- expression_matrix["Tnfaip3", ]
  
  correlations <- apply(expression_matrix, 1, function(gene_expression) {
    cor(tnfaip3_expression, gene_expression, use = "pairwise.complete.obs")
  })
  
  correlation_df <- data.frame(Gene = rownames(expression_matrix), CorrelationWithTnfaip3 = correlations)
  correlation_df_sorted <- correlation_df[order(-abs(correlation_df$CorrelationWithTnfaip3)), ]
  
  return(correlation_df_sorted)
}

# Calculate correlations for each condition
correlations_Lean <- calculate_correlations(Lean)
correlations_Obese <- calculate_correlations(Obese)
correlations_WL <- calculate_correlations(WL)

# Adding condition column for identification
correlations_Lean$Condition <- 'Lean'
correlations_Obese$Condition <- 'Obese'
correlations_WL$Condition <- 'WL'

# Combine all correlations into one dataframe
all_correlations <- rbind(correlations_Lean, correlations_Obese, correlations_WL)
library(ComplexHeatmap)

# Assuming you want to plot correlations of top N genes with highest variance in correlations across conditions
# This step is optional and can be customized
top_genes <- head(sort(apply(all_correlations[,-(1:2)], 1, var), decreasing = TRUE), 20)

heatmap_data <- all_correlations[all_correlations$Gene %in% top_genes, ]

# Create the heatmap
Heatmap(as.matrix(t(heatmap_data[,-1])),
        name = "Correlation",
        row_labels = heatmap_data$Gene,
        column_title = "Correlation with TNFAIP3 by Condition",
        cluster_columns = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE)

# Convert to data frame if it's not already one
all_correlations <- data.frame(all_correlations)

# Ensure CorrelationWithTnfaip3 is numeric for sorting
all_correlations$CorrelationWithTnfaip3 <- as.numeric(as.character(all_correlations$CorrelationWithTnfaip3))

# Sorting data by Condition and CorrelationWithTnfaip3
all_correlations <- all_correlations[order(all_correlations$Condition, -abs(all_correlations$CorrelationWithTnfaip3)),]

library(dplyr)

# Function to select top 5 positive and negative correlations
select_top_genes <- function(data) {
  top_positive <- data %>% 
    filter(CorrelationWithTnfaip3 > 0) %>% 
    top_n(5, CorrelationWithTnfaip3)
  
  top_negative <- data %>% 
    filter(CorrelationWithTnfaip3 < 0) %>% 
    top_n(-5, CorrelationWithTnfaip3)
  
  bind_rows(top_positive, top_negative)
}

# Applying the function to each condition and combining results
top_genes_by_condition <- all_correlations %>% 
  group_by(Condition) %>% 
  do(select_top_genes(.))

# Removing duplicate genes across conditions
top_unique_genes <- top_genes_by_condition %>% 
  distinct(Gene, .keep_all = TRUE)

library(ggplot2)
library(ComplexHeatmap)

# Preparing data for the heatmap (ensure it's in a matrix format where genes are rows and conditions are columns)
heatmap_data <- matrix(data = NA, nrow = length(unique(top_unique_genes$Gene)), ncol = length(unique(top_unique_genes$Condition)),
                       dimnames = list(unique(top_unique_genes$Gene), unique(top_unique_genes$Condition)))

# Filling the matrix with correlation values
for (i in 1:nrow(top_unique_genes)) {
  gene <- top_unique_genes$Gene[i]
  condition <- top_unique_genes$Condition[i]
  correlation <- top_unique_genes$CorrelationWithTnfaip3[i]
  heatmap_data[gene, condition] <- correlation
}

# Creating the heatmap
Heatmap(heatmap_data, name = "Correlation", 
        row_title = "Gene", column_title = "Condition",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE)

# Assuming heatmap_data is your matrix prepared for the heatmap

# Remove rows and columns with only NA values (if applicable)
heatmap_data <- heatmap_data[rowSums(is.na(heatmap_data)) < ncol(heatmap_data),]
heatmap_data <- heatmap_data[, colSums(is.na(heatmap_data)) < nrow(heatmap_data)]

# Option 1: Impute NA/NaN values with 0
heatmap_data[is.na(heatmap_data) | is.nan(heatmap_data)] <- 0

# Option 2: Alternatively, if you prefer to remove rows/genes with any NA/NaN values
# heatmap_data <- heatmap_data[!rowSums(is.na(heatmap_data) | is.nan(heatmap_data)),]
library(ComplexHeatmap)

Heatmap(heatmap_data, name = "Correlation", 
        row_title = "Gene", column_title = "Condition",
        cluster_rows = TRUE, cluster_columns = TRUE,
        show_row_names = TRUE, show_column_names = TRUE,
        na_col = "white") # This option colors NA values as white (if any are left)

# Find markers and limit to those expressed in greater than 75% of target population
all_markers <- FindAllMarkers(object = df) %>%
  Add_Pct_Diff() %>%
  filter(pct_diff > 0.6)

top_markers <- Extract_Top_Markers(marker_dataframe = all_markers, num_genes = 7, named_vector = FALSE,
                                   make_unique = TRUE)

Clustered_DotPlot(seurat_object = df, features = top_markers)

#Back to corrs:
library(ggplot2)

# Assuming 'top_genes_by_condition' is your data frame
top_genes_by_condition_filtered <- top_genes_by_condition[top_genes_by_condition$Gene != 'Tnfaip3', ]

library(ggplot2)

library(ggplot2)

ggplot(top_genes_by_condition_filtered, aes(x = Gene, y = CorrelationWithTnfaip3, color = Condition)) +
  geom_point(size = 5) + # Increase dot size
  scale_color_manual(values = c("Lean" = "black", "Obese" = "red", "WL" = "steelblue")) + # Custom colors
  facet_wrap(~ Condition, scales = "free_x") +
  theme_minimal() +
  labs(title = "Gene Correlation with TNFAIP3 Across Diets",
       x = "Gene",
       y = "Correlation with TNFAIP3") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12), # Make X axis labels bold and increase size
    axis.text.y = element_text(face = "bold", size = 12), # Make Y axis labels bold and increase size
    strip.text = element_text(face = "bold"), # Make facet labels bold
    panel.grid.major = element_blank(), # Remove major gridlines
    panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.title = element_text(face = "bold", size = 14), # Bold and increase size of axis titles
    axis.ticks.x = element_line(color = "black") # Add vertical marks (ticks) on the x-axis
  ) +
  geom_vline(xintercept = seq(0.5, nrow(top_genes_by_condition_filtered) + 0.5, by = 1), linetype = "dotted", color = "grey") # Add vertical dotted lines for each gene


#correlations:
library(dplyr)
library(Matrix)

pseudobulk_by_celltype_manual <- function(seurat_object) {
  # Get the data matrix
  data_matrix <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
  
  # Ensure that the clusters column is correctly set as factor
  seurat_object$clusters <- as.factor(seurat_object$clusters)
  
  # Create a dataframe linking cell names to clusters
  cell_clusters <- data.frame(cell = names(seurat_object$clusters),
                              cluster = seurat_object$clusters,
                              stringsAsFactors = FALSE)
  
  # Initialize an empty list to store pseudobulked data
  pseudobulked_data <- list()
  
  # Loop through each cluster, aggregate expression data, and store in list
  for(cluster in unique(cell_clusters$cluster)) {
    cells_in_cluster <- cell_clusters$cell[cell_clusters$cluster == cluster]
    
    # Subset data matrix to cells in current cluster and calculate mean expression
    cluster_data <- data_matrix[, cells_in_cluster]
    mean_expression <- rowMeans(cluster_data)
    
    # Store mean expression in list
    pseudobulked_data[[cluster]] <- mean_expression
  }
  
  # Combine the pseudobulked data into a single matrix
  combined_data <- do.call(cbind, pseudobulked_data)
  
  # Convert to a Seurat object if needed, or proceed with combined_data as is
  return(combined_data)
}

# Apply manual pseudobulking
Lean_pseudobulk_data <- pseudobulk_by_celltype_manual(Lean)
Obese_pseudobulk_data <- pseudobulk_by_celltype_manual(Obese)
WL_pseudobulk_data <- pseudobulk_by_celltype_manual(WL)

calculate_correlations_pseudobulk <- function(pseudobulk_data) {
  # Ensure TNFAIP3 is in the rownames
  if(!"TNFAIP3" %in% rownames(pseudobulk_data)) {
    stop("TNFAIP3 not found in the dataset")
  }
  
  # Extract TNFAIP3 expression
  tnfaip3_expression <- pseudobulk_data["TNFAIP3", ]
  
  # Initialize a vector to store correlations
  correlations <- numeric(nrow(pseudobulk_data) - 1)
  
  # Loop through each gene and calculate its correlation with TNFAIP3
  for(i in seq_along(correlations)) {
    gene_expression <- pseudobulk_data[i, ]
    correlations[i] <- cor(tnfaip3_expression, gene_expression, use = "complete.obs")
  }
  
  # Create a data frame of results
  correlation_df <- data.frame(
    Gene = rownames(pseudobulk_data)[-which(rownames(pseudobulk_data) == "TNFAIP3")],
    CorrelationWithTNFAIP3 = correlations
  )
  
  return(correlation_df)
}

# Note: Ensure that 'Lean_pseudobulk_data', 'Obese_pseudobulk_data', and 'WL_pseudobulk_data' are correctly formatted matrices
correlations_Lean <- calculate_correlations_pseudobulk(Lean_pseudobulk_data)
correlations_Obese <- calculate_correlations_pseudobulk(Obese_pseudobulk_data)
correlations_WL <- calculate_correlations_pseudobulk(WL_pseudobulk_data)

# Add condition labels for identification
correlations_Lean$Condition <- 'Lean'
correlations_Obese$Condition <- 'Obese'
correlations_WL$Condition <- 'WL'

# Combine all correlations into one dataframe
all_correlations <- rbind(correlations_Lean, correlations_Obese, correlations_WL)

ggplot(top_genes_by_condition, aes(x = Gene, y = CorrelationWithTnfaip3, color = Condition)) +
  geom_point(size = 5) + # Increase dot size for visibility
  scale_color_manual(values = c("Lean" = "black", "Obese" = "red", "WL" = "steelblue")) + # Assign custom colors
  facet_wrap(~ Condition, scales = "free_x", ncol = 1) + # Separate plot for each condition with free x scales
  theme_minimal() + # Use a minimalistic theme
  labs(title = "Top Gene Correlations with TNFAIP3 Across Diets", x = "Gene", y = "Correlation with TNFAIP3") + # Set plot title and axis labels
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12), # Make X axis labels bold and angled for readability
    axis.text.y = element_text(face = "bold", size = 12), # Make Y axis labels bold
    strip.text = element_text(face = "bold"), # Bold facet labels for clarity
    panel.grid.major = element_blank(), # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(), # Remove minor gridlines
    axis.title = element_text(face = "bold", size = 14), # Bold and enlarge axis titles
    axis.ticks.x = element_line(color = "black") # Enhance x-axis ticks
  ) +
  geom_vline(xintercept = seq_along(unique(top_genes_by_condition$Gene)) - 0.5, linetype = "dotted", color = "grey") # Add vertical lines for gene alignment

# Assuming 'all_correlations' is already created as per your previous step
# Remove TNFAIP3
all_correlations_filtered <- all_correlations[all_correlations$Gene != "TNFAIP3", ]


# Filter all_correlations_filtered to keep only the genes in formatted_genes
selected_gene_correlations <- all_correlations_filtered %>%
  filter(Gene %in% formatted_genes)

library(ggplot2)

ggplot(selected_gene_correlations, aes(x = Gene, y = CorrelationWithTnfaip3, color = Condition)) +
  geom_point(size = 4) +  # Increase dot size for visibility
  scale_color_manual(values = c("Lean" = "black", "Obese" = "red", "WL" = "steelblue")) +  # Custom colors for each diet
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at y=0
  theme_minimal() +  # Use a minimalistic theme
  labs(title = "NFKB Gene Target Correlations with TNFAIP3 Across Diets", 
       x = "Gene", 
       y = "Correlation with TNFAIP3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),  # Bold and angled X axis labels for readability
        axis.text.y = element_text(face = "bold", size = 12),  # Bold Y axis labels
        strip.text = element_text(face = "bold"),  # Bold facet labels
        panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.title = element_text(face = "bold", size = 14))  # Bold and enlarged axis titles

# Assuming all_correlations_filtered is already created and formatted_genes is defined
# First, filter for the genes of interest
selected_gene_correlations <- all_correlations_filtered %>%
  filter(Gene %in% formatted_genes)

# Now, remove rows where correlation is between -0.1 and 0.1
selected_gene_correlations <- selected_gene_correlations %>%
  filter(CorrelationWithTnfaip3 <= -0.1 | CorrelationWithTnfaip3 >= 0.1)

# Proceed to plotting
ggplot(selected_gene_correlations, aes(x = Gene, y = CorrelationWithTnfaip3, color = Condition)) +
  geom_point(size = 5) +  # Increase dot size for visibility
  scale_color_manual(values = c("Lean" = "black", "Obese" = "red", "WL" = "steelblue")) +  # Custom colors for each diet
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at y=0
  theme_minimal() +  # Use a minimalistic theme
  labs(title = "NFKB Gene Target Correlations with TNFAIP3 Across Diets", 
       x = "Gene", 
       y = "Correlation with TNFAIP3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),  # Bold and angled X axis labels for readability
        axis.text.y = element_text(face = "bold", size = 12),  # Bold Y axis labels
        strip.text = element_text(face = "bold"),  # Bold facet labels
        panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.title = element_text(face = "bold", size = 14))  # Bold and enlarged axis titles




# Filter all_correlations_filtered to keep only the genes in formatted_genes
selected_gene_correlations <- all_correlations_filtered %>%
  filter(Gene %in% formatted_genes)

# Arrange the genes based on their expression in lean mice (assuming expression values are in a column called 'Expression_lean')
selected_gene_correlations <- selected_gene_correlations %>%
  arrange(desc(Expression_lean))  # Use 'desc' to sort in descending order

library(ggplot2)

ggplot(selected_gene_correlations, aes(x = Gene, y = CorrelationWithTnfaip3, color = Condition)) +
  geom_point(size = 4) +  # Increase dot size for visibility
  scale_color_manual(values = c("Lean" = "black", "Obese" = "red", "WL" = "steelblue")) +  # Custom colors for each diet
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Add horizontal line at y=0
  theme_minimal() +  # Use a minimalistic theme
  labs(title = "NFKB Gene Target Correlations with TNFAIP3 Across Diets", 
       x = "Gene", 
       y = "Correlation with TNFAIP3") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),  # Bold and angled X axis labels for readability
        axis.text.y = element_text(face = "bold", size = 12),  # Bold Y axis labels
        strip.text = element_text(face = "bold"),  # Bold facet labels
        panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        axis.title = element_text(face = "bold", size = 14))  # Bold and enlarged axis titles

library(monocle3)

#ZFP36 test
poscells <- WhichCells(df, expression = Zfp36l1 > 0.75)
df@meta.data$Zfp36_exp <- ifelse(colnames(df) %in% poscells, "Pos", "Neg")

gene_markers <- FindMarkers(df, ident.1 = 'Pos', ident.2 = 'Neg', group.by = 'Zfp36_exp', logfc.threshold = 0, max.cells.per.ident = 500)

genes_of_interest <- c("Nfkb1",    # NF-??B1 p50 subunit, part of the NF-??B complex
                       "Rela",     
                       "Ikbkb",    
                       "Ikbkg",    
                       "Tnf",  'Hif1a',  
                       "Il1b",    
                       "Ripk1",  
                       "Traf6",   
                       "Rel", 'Birc3', 'Map3k8', 'Ikbke',
                       "Tnfrsf1a",
                       "Tnfrsf1b", 
                       "Nfkbia", "Nfkb1", 'Relb', 'Nfkb2', 
                       'Chuk', 'Traf6', 'Traf5', 'Traf3', 
                       'Traf4', 'Il17ra', 'Il17rc', 'Tradd',
                       'Tax1bp1', 'Tab1', 'Tab2', 'Tab3', 'Tnip1', 'Tnip2', 'Tnip3', 'Il6', 'Nfkbiz', 'Nfkbid', 'Ikbke')  

gene_markers$X <- rownames(gene_markers)
#attention <- subset(gene_markers, X %in% genes_of_interest)
attention <- subset(gene_markers, X == genes_of_interest)

gene_markers$pvalue <- gene_markers$p_val_adj
gene_markers$log2FoldChange <- gene_markers$avg_log2FC
gene_markers$genes <- gene_markers$X

genes_to_remove <- c('Zfp36l1')

# Filter the data frame
gene_markers <- gene_markers[!gene_markers$X %in% genes_to_remove, ]
attention_genes <- gene_markers[gene_markers$X %in% genes_of_interest, ]

#attention <- subset(gene_markers, X %in% genes_of_interest)
#attention <- subset(gene_markers, X == genes_of_interest)

ggvolc(gene_markers, attention_genes, fc = 0.5, p_value = 0.05, add_seg = T) +
  labs(title="NF\u03BAB Genes - Zfp36 Positive vs Negative T Cells")

View(gene_markers)



######

markers <- FindMarkers(df, 
                       ident.1 = 'Lean', 
                       ident.2 = 'Obese',
                       group.by = 'Condition',
                       logfc.threshold = 0.25)

# View the markers
View(markers)

# Add gene names to the markers data frame
markers$genes <- rownames(markers)

# Subset for significant markers
pos_markers <- subset(markers, p_val_adj < 0.5 & avg_log2FC > 0.1)
neg_markers <- subset(markers, p_val_adj < 0.5 & avg_log2FC < -0.1)

# Extract gene names for pathway analysis
gene_names <- pos_markers$genes
print(gene_names)
# Perform pathway enrichment analysis
enrich_result <- enrichGO(gene = gene_names,
                          keyType = 'SYMBOL',
                          OrgDb = orgDb, pvalueCutoff = 0.05)
print(enrich_result)
# Visualize enriched pathways with a dot plot
dot_plot <- dotplot(enrich_result, showCategory = 20) +
  ggtitle(paste(ident1, 'vs', ident2, 'Enriched Pathways'))

# Return the dot plot and the enrichment result
list(dot_plot = dot_plot, enrich_result = enrich_result)

nfkb_receptor_genes_mouse <- c("Tnfrsf1a", # TNF receptor superfamily, member 1a (TNF alpha receptor)
                               "Tnfrsf1b", # TNF receptor superfamily, member 1b (TNF beta receptor)
                               "Cd40",     # CD40 molecule, TNF receptor superfamily member 5
                               "Ltbr",     # Lymphotoxin beta receptor
                               "Tnfrsf8",  # TNF receptor superfamily, member 8 (CD30)
                               "Tnfrsf9",  # TNF receptor superfamily, member 9 (4-1BB)
                               "Tnfrsf13c",# B cell-activating factor receptor
                               "Tnfrsf17", # TNF receptor superfamily, member 17 (BCMA)
                               "Tnfrsf13b",# TNF receptor superfamily, member 13B (TACI)
                               "Cd27",     # CD27 molecule
                               "Traf3",    # TNF receptor-associated factor 3
                               "Traf2",    # TNF receptor-associated factor 2
                               "Nfkb1",    # Nuclear factor kappa B subunit 1
                               "Nfkb2",    # Nuclear factor kappa B subunit 2
                               "Rela",     # RELA proto-oncogene, NF-KB subunit
                               "Relb")     # RELB proto-oncogene, NF-KB subunit
genes <- c("Acadvl", "Hadha", "Hadhb", "Fabp3", "Fabp4", "Fabp5", "Gpam", "Agpat2")
VlnPlot(df, features = genes, stack = T, flip = T)
# Print the vector
print(nfkb_receptor_genes_mouse)

cell_labels <- Idents(df)

# Convert to a named vector
cell_labels_vector <- as.vector(cell_labels)
names(cell_labels_vector) <- names(cell_labels)

# Add this information to the Seurat object metadata
df <- AddMetaData(df, metadata = cell_labels_vector, col.name = "cell_type")

# Check the metadata to confirm the column was added
View(df@meta.data)
