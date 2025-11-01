# EXTENDED DATA FIGURE 6D
path = 'data_repo/' 
dat = readRDS(paste0(path, 'merged_post.rds'))

meta = dat@meta.data
meta$celltype = ''

Tcell = WhichCells(dat, expression = CD3E > 0.5 & (TRBC1 > 0.5 | TRDC > 0.5))
NK = WhichCells(dat, expression = CD3E < 0.5 & NCR1 > 0.5 & (FCGR3A > 0.5 | NCAM1 > 0.5))
APC = WhichCells(dat, expression = CD3E < 0.5 & NCR1 < 0.5 & (CD14 > 0.5 | ITGAM > 0.5 | ITGAX > 0.5) & (CD80 > 0.5 | CD86 > 0.5 | 'HLA-DRB1' > 0.5))
Bcell = WhichCells(dat, expression = IGHM > 0.5 & (CXCR5 > 0.5 | CD74 > 0.5))

meta[Tcell,]$celltype = 'T cell'
meta[NK,]$celltype = 'NK'
meta[APC,]$celltype = 'APC'
meta[Bcell,]$celltype = 'Bcell'

dat@meta.data = meta
prop.table(table(dat$celltype))
prop.table(table(dat$celltype, dat$POD), margin = 2)

# EXTENDED DATA FIGURE 6E
prop.table(table(dat$celltype, dat$seurat_clusters2), margin = 2)

# EXTENDED DATA FIGURE 6F
# Copy to clipboard: 2033 x 850
FeaturePlot(dat, features = c('CD69', 'ICOS', 'IFNG', 'STAT4', 'TBX21', 'GATA3', 'CCL4', 'CCL5', 'TNF', 'CD28'), ncol = 5, pt.size = 2, order = T)  &
  labs(x = "UMAP_1", y = "UMAP_2") &
  theme(axis.title  = element_text(size = 24),
        axis.text   = element_text(size = 22),
        legend.title= element_text(size = 24),
        legend.text = element_text(size = 16),
        plot.title   = element_text(face = "italic", size = 24))

# EXTENDED DATA FIGURE 6G
FeaturePlot(dat, features = c('CD68', 'SIRPA', 'MRC1', 'FCGR3A', 'ITGAX', 'FLT3', 'TLR4'), ncol = 5) &
  labs(x = "UMAP_1", y = "UMAP_2") &
  theme(axis.title  = element_text(size = 24),
        axis.text   = element_text(size = 22),
        legend.title= element_text(size = 24),
        legend.text = element_text(size = 16),
        plot.title   = element_text(face = "italic", size = 24))
