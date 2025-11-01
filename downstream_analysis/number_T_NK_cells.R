path = 'data_repo/' 
dat = readRDS(paste0(path, 'merged_post.rds'))

Idents(dat) = 'POD'

# POD 33 T and NK cell count
meta = dat@meta.data
Tcell = WhichCells(dat, expression = CD3E > 0.5 & (TRBC1 > 0.5 | TRDC > 0.5), idents = 'POD33K') # 92
NK = WhichCells(dat, expression = CD3E < 0.5 & NCR1 > 0.5 & (FCGR3A > 0.5 | NCAM1 > 0.5), idents = 'POD33K') # 48

# TCRb sequences by POD
meta = dat@meta.data
meta_POD14K = meta %>% filter(POD == 'POD14K')
meta_POD14L = meta %>% filter(POD == 'POD14L')
meta_POD28K = meta %>% filter(POD == 'POD28K')
meta_POD33K = meta %>% filter(POD == 'POD33K')
meta_POD61K = meta %>% filter(POD == 'POD61K')

res = sum(meta_POD14K$TCRb_1 != '0') # 0
res = sum(meta_POD14L$TCRb_1 != '0') # 1341
res = sum(meta_POD28K$TCRb_1 != '0') # 7
res = sum(meta_POD33K$TCRb_1 != '0') # 73
res = sum(meta_POD61K$TCRb_1 != '0') # 8
