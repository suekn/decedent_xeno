# EXTENDED DATA FIGURE 6A

# Hybrid Reference
library(Seurat)

path = 'data_repo/' 
POD14k = Read10X(data.dir = paste0(path, 'KB03'))
POD14L = Read10X(data.dir = paste0(path, 'KB02'))
POD28 = Read10X(data.dir = paste0(path, 'KB04'))
POD33 = Read10X(data.dir = paste0(path, 'KB05'))
POD61 = Read10X(data.dir = paste0(path, 'KB06'))

POD14k = CreateSeuratObject(counts = POD14k, min.features = 500, min.cells = 15)
POD14L = CreateSeuratObject(counts = POD14L, min.features = 500, min.cells = 15)
POD28 = CreateSeuratObject(counts = POD28, min.features = 500, min.cells = 10)
POD33 = CreateSeuratObject(counts = POD33, min.features = 500, min.cells = 15)
POD61 = CreateSeuratObject(counts = POD61, min.features = 300, min.cells = 15)

process_dataset <- function(dataset, pattern = '^MT-', percent_mt_col = 'percent.mt', min_ncount = 1000, max_ncount = 25000) {
  dataset <- PercentageFeatureSet(dataset, pattern = pattern, col.name = percent_mt_col)
  dataset <- subset(dataset, subset = dataset[[percent_mt_col]] < 15 & nCount_RNA > min_ncount & nCount_RNA < max_ncount)
  return(dataset)
}

POD14k <- process_dataset(POD14k)
POD14L <- process_dataset(POD14L)
POD28 <- process_dataset(POD28)
POD33 <- process_dataset(POD33, max_ncount = 30000)
POD61 <- process_dataset(POD61, max_ncount = 30000)

POD14k@meta.data$POD = 'POD14K'
POD14L@meta.data$POD = 'POD14L'
POD28@meta.data$POD = 'POD28K'
POD33@meta.data$POD = 'POD33K'
POD61@meta.data$POD = 'POD61K'

POD14k@meta.data$rejection = 'non-rejecting'
POD14L@meta.data$rejection = 'non-rejecting'
POD28@meta.data$rejection = 'non-rejecting'
POD33@meta.data$rejection  = 'rejecting'
POD61@meta.data$rejection  = 'non-rejecting'

merged <- merge(POD14k, y = list(POD14L, POD28, POD33, POD61))
DefaultAssay(merged) = 'RNA'
merged = JoinLayers(merged)

RNA = merged@assays$RNA$counts %>% as.data.frame()
human_genes <- grep("^GRCh38-cellrangerV7", rownames(RNA), value = TRUE)
pig_genes <- grep("^susscrofa", rownames(RNA), value = TRUE)

human_counts <- RNA[human_genes, ]
pig_counts <- RNA[pig_genes, ]

human_counts_sum <- colSums(human_counts)
pig_counts_sum <- colSums(pig_counts)

summary_counts <- data.frame(
  Cell = colnames(RNA),
  Human_Counts = human_counts_sum,
  Pig_Counts = pig_counts_sum
)


library(ggplot2)
summary_counts$Above_Line <- ifelse(summary_counts$Pig_Counts > summary_counts$Human_Counts, "Above", "Below")

pig_cells = summary_counts %>% filter(Above_Line == 'Above') %>% row.names()
human_cells = summary_counts %>% filter(Above_Line == 'Below') %>% row.names()

above_below_counts <- table(summary_counts$Above_Line)
num_above <- above_below_counts["Above"]
num_below <- above_below_counts["Below"]

above_label <- paste("Above (Pig-dominant):", num_above)
below_label <- paste("Below (Human-dominant):", num_below)

# Create the scatterplot with the y = x line and add counts as annotations
# Copy from clipboard: 750 x 600 
ggplot(summary_counts, aes(x = Human_Counts, y = Pig_Counts)) +
  geom_point(aes(color = Above_Line), alpha = 0.6) +  # Color points based on classification
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 1) +  # Add y = x line
  theme_minimal() +                                # Use a clean theme
  labs(
    title = "Pig vs Human Cells",
    x = "Number of Human Gene Counts",
    y = "Number of Pig Gene Counts",
    color = "Position"
  ) +
    theme(
    plot.title  = element_text(size = 30, hjust = 0.5),
    axis.title  = element_text(size = 26),
    axis.text   = element_text(size = 22),
    legend.title= element_text(size = 24),
    legend.text = element_text(size = 20)
  )

#write.csv(pig_cells, file = 'pig_cells.csv')
#write.csv(human_cells, file = 'human_cells.csv')

# Human Only Reference
dat = readRDS('merged_pre.rds')
meta = dat@meta.data
DimPlot(dat, group.by = 'seurat_clusters2')

pig_cells = read.csv('pig_cells.csv')

remove = intersect(row.names(meta), pig_cells$x)

sub = subset(dat, cells = remove, invert = T)

unique(sub$seurat_clusters2) # cluster 5 is removed
Idents(sub) <- sub$seurat_clusters2
sub <- RenameIdents(sub, `6` = "5", `7` = "6", `8` = "7", `9` = "8", `10` = "9", `11` = "10", `12` = "11")
Idents(sub) <- factor(as.character(Idents(sub)), levels = as.character(0:11)) # relevel clusters from 0 to 11
sub$seurat_clusters2 <- Idents(sub) 


# EXTENDED DATA FIGURE 6B
# Copy from clipboard: 650 x 600 
xlim <- range(-14:22)
ylim <- range(-14:22)

DimPlot(dat, group.by = 'seurat_clusters2', label = F, pt.size = 2) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))


# EXTENDED DATA FIGURE 6C
DimPlot(sub, group.by = 'seurat_clusters2', label = F, pt.size = 2) + NoLegend() + ggtitle(NULL) + coord_equal(xlim = xlim, ylim = ylim) + labs(x = "UMAP_1", y = "UMAP_2") +  theme(axis.title = element_text(size = 30), 
                   axis.text  = element_text(size = 28),
                   legend.title = element_text(size = 30),
                   legend.text  = element_text(size = 28))


saveRDS(sub, file = 'merged_post.rds')
