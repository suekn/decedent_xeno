dat = readRDS('merged_post.rds')

# FIGURE 4A
Idents(dat) = 'POD'

xlim <- range(-14:22)
ylim <- range(-14:22)

POD14K = subset(dat, idents = "POD14K")
POD14L = subset(dat, idents = "POD14L")
POD28K = subset(dat, idents = "POD28K")
POD33K = subset(dat, idents = "POD33K")
POD61K = subset(dat, idents = "POD61K")

# Copy from clipboard: 500 x 450 
colors <- setNames(c("#F8766D", "#DE8C00", "#B79F00", "#7CAE00", "#1AC14C", "#00C08B", "#00BFC4", "#0DB8F1", "#619CFF", "#C77CFF", "#F564E3", "#FF64B0"), as.character(0:11))

DimPlot(dat, group.by = 'seurat_clusters2', pt.size = 2, cols = colors) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))

DimPlot(POD14K, group.by = 'seurat_clusters2', pt.size = 2, cols = colors) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))

DimPlot(POD14L, group.by = 'seurat_clusters2', pt.size = 2, cols = colors) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))

DimPlot(POD28K, group.by = 'seurat_clusters2', pt.size = 2, cols = colors) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))

DimPlot(POD33K, group.by = 'seurat_clusters2', pt.size = 2, cols = colors) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))

meta = dat@meta.data
prop.table(table(meta$seurat_clusters2))
prop.table(table(meta$seurat_clusters2, meta$POD), margin = 2)


# FIGURE 4B
meta = dat@meta.data
meta$celltype = ''

Tcell = WhichCells(dat, expression = CD3E > 0.5 & (TRBC1 > 0.5 | TRDC > 0.5))
NK = WhichCells(dat, expression = CD3E < 0.5 & NCR1 > 0.5 & (FCGR3A > 0.5 | NCAM1 > 0.5))
APC = WhichCells(dat, expression = CD3E < 0.5 & NCR1 < 0.5 & (CD14 > 0.5 | ITGAM > 0.5 | ITGAX > 0.5) & (CD80 > 0.5 | CD86 > 0.5 | 'HLA-DRB1' > 0.5))
Bcell = WhichCells(dat, expression = IGHM > 0.5 & (CXCR5 > 0.5 | CD74 > 0.5))

meta[Tcell,]$celltype = 'T cell'
meta[NK,]$celltype = 'NK'
meta[APC,]$celltype = 'APC'
meta[Bcell,]$celltype = 'B cell'

dat@meta.data = meta
colors = setNames(c("grey80", "#00BFC4", "#F8766D", "#7CAE00", "#DE8C00"), c("", "T cell", "APC", "NK", "B cell"))
layers = c("NK", "T cell", "APC", "B cell", "")
DimPlot(dat, group.by = 'celltype', pt.size = 1.5, cols = colors, order = layers) + NoLegend() + ggtitle(NULL) + labs(x = "UMAP_1", y = "UMAP_2") +  
  theme(axis.title = element_text(size = 30), 
        axis.text  = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 28))


# FIGURE 4C
library(ComplexHeatmap)
library(circlize)

genes = c('CD2', 'CD3D', 'CD3E', 'CD3G', 'CD4', 'CD8A', 'TRDC', 'CD28', 'TBX21', 'GATA3', 'FOXP3', 'BCL6', 'ASCL2', 'ICOS', 'GNLY', 'CCL5', 'CCL4', 'CXCL8', 'CXCL11', 'GZMB', 'PRF1', 'IL2RA', 'TNF', 'CCR7', 'CD69', 'CD274', 'FCGR1A', 'FCGR1B', 'FCGR2A', 'FCGR2B', 'FCGR3A', 'NCR1', 'NCR3', 'KLRC1', 'KLRC2', 'NCAM1', 'KLRK1', 'CD14', 'CD80', 'CD86', 'ITGAX', 'ITGAM', 'CD163', 'IFNG', 'HLA-A', 'HLA-B', 'HLA-DRB1')

# Assuming heatmap_data, annotation_col, and selected_genes are already defined
# Create the color palette
selected_genes = genes
dat = ScaleData(dat, features = rownames(dat))
heatmap_data <- as.matrix(GetAssayData(dat, slot = "scale.data")[selected_genes, ])

Idents(dat) = dat$seurat_clusters2
clusters <- Idents(dat)
col_order <- order(clusters)
heatmap_data <- heatmap_data[, col_order] # order heatmap by clusters

# Assuming heatmap_data and col_order are defined
clusters <- Idents(dat)[col_order]
rejection_status <- dat$rejection[col_order]

# Combine cluster and rejection status into a data frame
metadata <- data.frame(Cluster = clusters, Rejection = rejection_status)

# Order by Rejection, then by Cluster
ordered_metadata <- metadata %>%
  arrange(Rejection, Cluster)

# Reorder heatmap_data according to the new order
heatmap_data <- heatmap_data[, rownames(ordered_metadata)]
clusters <- ordered_metadata$Cluster
rejection_status <- ordered_metadata$Rejection

# Create a palette for clusters
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
color_list <- ggplotColours(n=12)

colors <- color_list
cluster_colors <- setNames(colors, levels(as.factor(clusters)))

colors <- c("cornflowerblue", "coral3")
rejection_colors = setNames(colors, levels(as.factor(rejection_status)))
ha <- HeatmapAnnotation(Cluster = as.factor(clusters), 
                        Rejection = factor(rejection_status),       
                        col = list(Cluster = cluster_colors, Rejection = rejection_colors),
                        annotation_name_side = 'left')

x <- as.numeric(heatmap_data)
lims <- quantile(x, c(.02, .98), na.rm = TRUE) # avoid extremes to avoid washing out data
M <- max(abs(lims)) # any values greater than limits will be clipped at top color

col_fun <- colorRamp2(
  c(-M, -M*0.5, 0, M*0.5, M),
  c("#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#B2182B")
)

# Copy to clipboard: 700 x 907
Heatmap(heatmap_data, col = col_fun, cluster_columns = FALSE, cluster_rows = FALSE, row_names_gp = gpar(fontface = "italic"),
              show_column_names = FALSE, top_annotation = ha, show_row_names = T)



# Figure 4D
Idents(dat) = 'POD'

xlim <- range(-14:22)
ylim <- range(-14:22)

POD14K = subset(dat, idents = "POD14K")
POD14L = subset(dat, idents = "POD14L")
POD28K = subset(dat, idents = "POD28K")
POD33K = subset(dat, idents = "POD33K")
POD61K = subset(dat, idents = "POD61K")

# Copy from clipboard: 500 x 450 
colors = setNames(c("grey80", "red", "blue"), c("", "XDRTCC", "nonXDRTCC"))
layers = c("XDRTCC", "nonXDRTCC", "")

DimPlot(POD14K, group.by = 'XDRTCC', pt.size = 2, cols = colors, order = layers) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  
  theme(axis.title = element_text(size = 30), 
        axis.text  = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 28))

DimPlot(POD14L, group.by = 'XDRTCC', pt.size = 2, cols = colors, order = layers) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  
  theme(axis.title = element_text(size = 30), 
        axis.text  = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 28))

DimPlot(POD28K, group.by = 'XDRTCC', pt.size = 2, cols = colors, order = layers) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  
  theme(axis.title = element_text(size = 30), 
        axis.text  = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 28))

DimPlot(POD33K, group.by = 'XDRTCC', pt.size = 2, cols = colors, order = layers) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  
  theme(axis.title = element_text(size = 30), 
        axis.text  = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 28))

DimPlot(POD61K, group.by = 'XDRTCC', pt.size = 2, cols = colors, order = layers) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  
  theme(axis.title = element_text(size = 30), 
        axis.text  = element_text(size = 28),
        legend.title = element_text(size = 30),
        legend.text  = element_text(size = 28))


# Figure 4E
Idents(dat) = "seurat_clusters2"
dat = ScaleData(dat, features = rownames(dat@assays$RNA))

# marker selection
markers = FindAllMarkers(dat)
all = markers %>% group_by(cluster) %>% top_n(70, wt = avg_log2FC) %>% filter(cluster %in% c("3","5"))
top    = all %>% filter(cluster == "3") %>% top_n(28, wt = avg_log2FC)
bottom = all %>% filter(cluster == "5") %>% top_n(32, wt = avg_log2FC)

select = c(top$gene, bottom$gene)
sub = subset(dat, idents = c("3","5"))

# matrix in the order of 'select'
genes = select[!duplicated(select)]
genes = genes[genes %in% rownames(sub)]
mat   = GetAssayData(sub, assay = DefaultAssay(sub), slot = "scale.data")[genes, , drop = FALSE]

# clipping limits and colors
x    = as.numeric(mat)
lims = stats::quantile(x, c(.02, .98), na.rm = TRUE)
M    = max(abs(lims))
col_fun = colorRamp2(c(-M, 0, M), c("#2166AC", "#F7F7F7", "#B2182B"))

# order and annotate by cluster
clusters <- factor(Idents(sub), levels = c("3","5"))
ord <- order(clusters)
mat <- mat[, ord, drop = FALSE]                 # reorder matrix
clusters_ord <- clusters[ord]                   # reordered labels
cluster_col <- c("3" = "#CECBD0", "5" = "#BD98A2")
ha <- HeatmapAnnotation(Cluster = clusters_ord, col = list(Cluster = cluster_col), simple_anno_size = unit(3, "mm"), annotation_name_side = "left")

# Copy from clipboard: 520 x 907
Heatmap(mat, name = "Expression", col = col_fun, cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, show_row_names = TRUE, row_names_gp = gpar(col = "black", fontface = "italic"), column_split = clusters_ord, cluster_column_slices = FALSE, top_annotation = ha, use_raster = TRUE)



# Figure 4F
markers = FindMarkers(dat, ident.1 = '3', ident.2 = '5')
markers$gene = row.names(markers)

labels = c('NCAM', 'NCR1', 'NCR3', 'TRDC', 'TRGC1', 'TRGV9', 'TRGV7', 'GZMB', 'GZMK', 'GZMM', 'GZMA', 'GZMH', 'KLRK1', 'CD8B', 'CD3', 'CD3D', 'CD3E', 'CD3G', 'TRBC2')

# Copy from clipboard: 600 x 907
EnhancedVolcano(markers, lab = row.names(markers), x = 'avg_log2FC', y = 'p_val_adj', drawConnectors = TRUE, widthConnectors = 1.5, lengthConnectors = unit(0, "mm"), colConnectors = "black", pCutoff = 1e-5, labSize = 7, axisLabSize = 24, max.overlaps = 30, selectLab = labels) + ggtitle('Cluster 3 vs 5') 


# FIGURE 4G, H
source("DEEnrichr_NS2.R")
library(enrichR)

# P val adjust uses Benjamini Hochberg "Enrichr relies on several ranking methods to compute enrichment, which can be chosen by the user. The p-value is calculated with Fisher’s exact test, where genes are considered independent. The q-value is an adjusted p-value calculated using the Benjamini-Hochberg method." PMID: 33780170
# Copy from clipboard: 550 x 450 
DEenrichRPlot_NS(dat, ident.1 = '3', ident.2 = '5', max.genes = 65, balanced = F, enrich.database = "KEGG_2021_Human") 

# custom function uses FindMarkers 3 vs 5 (not FindAllMarkers)
DEenrichRPlot_NS(dat, ident.1 = '5', ident.2 = '3', max.genes = 65, balanced = F, enrich.database = "KEGG_2021_Human") 
