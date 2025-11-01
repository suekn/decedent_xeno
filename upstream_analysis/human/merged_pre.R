library(Seurat)
library(tidyverse)
library(SeuratObject)
library(cluster)
library(ggplot2)
library(EnhancedVolcano)

path = 'data_repo/' 
POD14k = Read10X(data.dir = paste0(path, 'KB03/filtered_feature_bc_matrix/'))
POD14L = Read10X(data.dir = paste0(path, 'KB02/filtered_feature_bc_matrix/' ))
POD28 = Read10X(data.dir = paste0(path, 'KB04/filtered_feature_bc_matrix/'))
POD33 = Read10X(data.dir = paste0(path, 'KB05/filtered_feature_bc_matrix/' ))
POD61 = Read10X(data.dir = paste0(path, 'KB06/filtered_feature_bc_matrix/' ))

POD14k = CreateSeuratObject(counts = POD14k, min.features = 500, min.cells = 50)
POD14L = CreateSeuratObject(counts = POD14L, min.features = 500, min.cells = 50)
POD28 = CreateSeuratObject(counts = POD28, min.features = 500, min.cells = 50)
POD33 = CreateSeuratObject(counts = POD33, min.features = 500, min.cells = 50)
POD61 = CreateSeuratObject(counts = POD61, min.features = 300, min.cells = 50)

POD14k = PercentageFeatureSet(POD14k, pattern = '^MT-', col.name = 'percent.mt')
POD14k = subset(POD14k, subset = percent.mt<15 & nCount_RNA>1000 & nCount_RNA <25000)

POD14L = PercentageFeatureSet(POD14L, pattern = '^MT-', col.name = 'percent.mt')
POD14L = subset(POD14L, subset = percent.mt<15 & nCount_RNA>1000 & nCount_RNA <25000)

POD28 = PercentageFeatureSet(POD28, pattern = '^MT-', col.name = 'percent.mt')
POD28 = subset(POD28, subset = percent.mt<15 & nCount_RNA>1000 & nCount_RNA <25000)

POD33 = PercentageFeatureSet(POD33, pattern = '^MT-', col.name = 'percent.mt')
POD33 = subset(POD33, subset = percent.mt<15 & nCount_RNA>1000 & nCount_RNA <30000)

POD61 = PercentageFeatureSet(POD61, pattern = '^MT-', col.name = 'percent.mt')
POD61 = subset(POD61, subset = percent.mt<15 & nCount_RNA>1000 & nCount_RNA <30000)

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
merged = NormalizeData(merged) 
merged = JoinLayers(merged)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged <- CellCycleScoring(merged, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

merged <- ScaleData(merged, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(merged))
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged = RunPCA(merged, features = VariableFeatures(object = merged))
merged = RunUMAP(merged, dims = 1:30, verbose = FALSE, umap.method = 'umap-learn', metric = 'correlation')
merged = FindNeighbors(merged, dims = 1:30, verbose = F)
merged = FindClusters(merged, resolution = seq(0.01,2, by= 0.1), verbose = F, algorithm = 1)
DimPlot(merged, label = T, group.by = 'RNA_snn_res.0.41')

merged$seurat_clusters2 = merged$RNA_snn_res.0.41
Idents(merged) = merged$seurat_clusters2
dat = merged

# Link TCR data
path = 'data_repo/' 
TCR_POD14L = read.csv(paste0(path, 'KB02T/filtered_contig_annotations.csv'))
TCR_POD28 = read.csv(paste0(path, 'KB04T/filtered_contig_annotations.csv'))
TCR_POD33 = read.csv(paste0(path, 'KB05T/filtered_contig_annotations.csv'))
TCR_POD61 = read.csv(paste0(path, 'KB06T/filtered_contig_annotations.csv'))

TCR_POD14L$pt = 'POD14L'
TCR_POD28$pt = 'POD28K'
TCR_POD33$pt = 'POD33K'
TCR_POD61$pt = 'POD61K'

TCR_POD14L$unique_barcode = paste0('POD14L', '_', TCR_POD14L$barcode)
TCR_POD28$unique_barcode = paste0('POD28K', '_', TCR_POD28$barcode)
TCR_POD33$unique_barcode = paste0('POD33K', '_', TCR_POD33$barcode)
TCR_POD61$unique_barcode = paste0('POD61K', '_', TCR_POD61$barcode)

# Check if there are identical barcodes between the different samples (this would throw off matching based on cell barcode in next step)
temp = sum(TCR_POD14L$barcode %in% TCR_POD28$barcode) #0
temp = sum(TCR_POD14L$barcode %in% TCR_POD33$barcode) #0
temp = sum(TCR_POD14L$barcode %in% TCR_POD61$barcode) #0
temp = sum(TCR_POD28$barcode %in% TCR_POD33$barcode) #0
temp = sum(TCR_POD14L$barcode %in% TCR_POD61$barcode) #0
temp = sum(TCR_POD33$barcode %in% TCR_POD61$barcode) #0

#tcr_seq = rbind(TCR_POD14L, TCR_POD28, TCR_POD33, TCR_POD61)
tcr_seq = rbind(TCR_POD14L, TCR_POD28, TCR_POD33, TCR_POD61)

dat@meta.data$barcodes = rownames(dat@meta.data)
dat@meta.data$barcodes = gsub("_\\d$", "", dat@meta.data$barcodes)
dat@meta.data$unique_barcode = paste0(dat@meta.data$POD, '_', dat@meta.data$barcodes) # generate unique barcodes by sample and barcode

# Assign TCRaa CDR3 sequence in metadata
dat@meta.data$TCRa_1 = 0
dat@meta.data$TCRa_2 = 0
dat@meta.data$TCRa_3 = 0
dat@meta.data$TCRb_1 = 0
dat@meta.data$TCRb_2 = 0
dat@meta.data$TCRb_3 = 0
l1 = nrow(dat@meta.data)

for(i in 1:l1){
  ind = which(tcr_seq$unique_barcode == dat@meta.data$unique_barcode[i])
  
  ind_2 = which(tcr_seq[ind,]$chain == 'TRA')
  l2 = length(ind_2)
  
  if(l2 == 1)
      dat@meta.data$TCRa_1[i] = tcr_seq[ind[ind_2[1]],]$cdr3
  if(l2 == 2){
      dat@meta.data$TCRa_1[i] = tcr_seq[ind[ind_2[1]],]$cdr3
      dat@meta.data$TCRa_2[i] = tcr_seq[ind[ind_2[2]],]$cdr3
    }
  if(l2 == 3){
      dat@meta.data$TCRa_1[i] = tcr_seq[ind[ind_2[1]],]$cdr3
      dat@meta.data$TCRa_2[i] = tcr_seq[ind[ind_2[2]],]$cdr3
      dat@meta.data$TCRa_3[i] = tcr_seq[ind[ind_2[3]],]$cdr3
    }
  
  
  ind_3 = which(tcr_seq[ind,]$chain == 'TRB')
  l3 = length(ind_3)
    if(l3 == 1)
      dat@meta.data$TCRb_1[i] = tcr_seq[ind[ind_3[1]],]$cdr3
    if(l3 == 2){
      dat@meta.data$TCRb_1[i] = tcr_seq[ind[ind_3[1]],]$cdr3
      dat@meta.data$TCRb_2[i] = tcr_seq[ind[ind_3[2]],]$cdr3
    }
    if(l3 == 3){
      dat@meta.data$TCRb_1[i] = tcr_seq[ind[ind_3[1]],]$cdr3
      dat@meta.data$TCRb_2[i] = tcr_seq[ind[ind_3[2]],]$cdr3
      dat@meta.data$TCRb_3[i] = tcr_seq[ind[ind_3[3]],]$cdr3
    }
}

# Identify XRTCCs
path = 'data_repo/'
DRTCC_CD8 = read.csv(paste0(path, 'MasterCD8DRTCCs_AgainstPreTx.csv'))
DRTCC_CD4 = read.csv(paste0(path, 'MasterCD4DRTCCs_AgainstPreTx.csv'))
nonDRTCC = read.csv(paste0(path, 'PreTxNonDRTCCsConservative.csv'))

meta = dat@meta.data

# XDRTCC
ind = which(meta$TCRb_1 %in% DRTCC_CD8$Amino.Acid | meta$TCRb_2 %in% DRTCC_CD8$Amino.Acid | meta$TCRb_3 %in% DRTCC_CD8$Amino.Acid) 
barcode = meta[ind,] %>% row.names()

ind = which(meta$TCRb_1 %in% DRTCC_CD4$Amino.Acid | meta$TCRb_2 %in% DRTCC_CD4$Amino.Acid | meta$TCRb_3 %in% DRTCC_CD4$Amino.Acid)
barcode2 = meta[ind,] %>% row.names()

highlight = c(barcode, barcode2)

dat$XDRTCC = ''
ind = which(row.names(dat@meta.data) %in% highlight)
dat$XDRTCC[ind] = 'XDRTCC'

# Non DRTCC
ind = which(meta$TCRb_1 %in% nonDRTCC$Amino.Acid | meta$TCRb_2 %in% nonDRTCC$Amino.Acid | meta$TCRb_3 %in% nonDRTCC$Amino.Acid)
barcode3 = meta[ind,] %>% row.names()

highlight2 = c(barcode3)

ind = which(row.names(dat@meta.data) %in% highlight2)
dat$XDRTCC[ind] = 'nonXDRTCC'

# XDRTCC CDR3
dat$XRTCC_cdr3 = ''
ind = which(dat$XDRTCC == 'XDRTCC')
dat$XRTCC_cdr3[ind] = dat$TCRb_1[ind]

saveRDS(dat, file = 'merged_pre.rds')
