###Pre-process
cancer <- subset(npc, subset = MajorCluster == "Malignant cells")
cancer_filt <- subset(cancer, subset = Source == "Tumor")
aa <- sample(colnames(cancer_filt),500)
cancer_filt_sample <- subset(cancer_filt, cells=aa)

###
cd8_filt_t <- subset(cd8_filt, subset = Source == "Tumor")
cd8_filt_sample <- subset(cd8_filt_t, cells=SeuratObject::WhichCells(cd8_filt_t, downsample=300))
cd8_filt_sample$annotation <- cd8_filt_sample@active.ident

###
subset_cells <-merge(x = cancer_filt_sample,y = list(cd8_filt_sample))

###
eff <- read.table("genes.results",
                  header = T,stringsAsFactors = F)
eff$gene <- strsplit(eff$gene_id,split = "_")  %>% lapply(.,function(x){x[[2]]})  %>% unlist()
eff_length <- plyr::mapvalues(rownames(subset_cells),from =eff$gene,to =eff$effective_length  )

ref_dat <- data.frame(gene = rownames(subset_cells), eff_length = eff_length,stringsAsFactors = F)
ref_dat$eff_length <- as.numeric(ref_dat$eff_length)
ref_dat <- na.omit(ref_dat)
ref_dat_rm <- ref_dat[ref_dat$eff_length > 0,]

###
dat <- subset_cells@assays$RNA@counts %>% as.matrix()
dat_filt <- dat[ref_dat_rm$gene,]
cpm <- RelativeCounts(data = dat_filt,scale.factor = 1e6)

########
TPMs <- apply(cpm,2,function(x){
  x/(ref_dat_rm$eff_length)
})


###
TPM <- round(TPMs,2)
labelData <- data.frame(cells = rownames(subset_cells@meta.data),labels = subset_cells@meta.data$annotation)  
save(TPM,labelData,file = "npc_cancer_CD8_sample_filt.RData")

###CSOmap
load("npc_cancer_CD8_sample_filt.RData")
library(CSOmapR)
library(CSOmapR.demo)
LR <- read.table("LR_pairs.txt",stringsAsFactors = F)

affinityMat = getAffinityMat(TPM, LR, verbose = T)

coords_res = runExactTSNE_R(
  X = affinityMat,
  no_dims = 3,
  max_iter = 1000,
  verbose = T
)
coords = coords_res$Y
rownames(coords) <- colnames(TPM)
colnames(coords) <- c('x', 'y', 'z')


require(dplyr)
# arrange data
coords_tbl = bind_cols(cellName = rownames(coords), as.data.frame(coords))

join_vec = setNames(colnames(labelData)[1], nm = colnames(coords_tbl)[1])
cellinfo_tbl = left_join(coords_tbl, labelData, by = join_vec)

density_obj = getDensity3D(cellinfo_tbl$x, cellinfo_tbl$y, cellinfo_tbl$z)
cellinfo_tbl = cellinfo_tbl %>% mutate(density = density_obj)

signif_results = getSignificance(coords, labels = cellinfo_tbl$labels, verbose = T)
contribution_list = getContribution(TPM, LR, signif_results$detailed_connections)

save(signif_results,contribution_list,cellinfo_tbl,file = "Cancer_CD8_res.RData")
