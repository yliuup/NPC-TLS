
NPC <- readRDS(paste0("data.rds"))
NPC <- SCTransform(NPC, assay = "Spatial", verbose = FALSE)
NPC <- RunPCA(NPC, assay = "SCT", verbose = FALSE)
NPC <- FindNeighbors(NPC, reduction = "pca", dims = 1:25)
NPC <- FindClusters(NPC, verbose = FALSE,resolution = 1,  random.seed = 88)
NPC <- RunUMAP(NPC, reduction = "pca", dims = 1:25, umap.method = "uwot",
              n.neighbors = 25L, min.dist = 0.5)

library(ggsci)
cols <- c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_lancet()(9))[-8]
png(paste0("NPC","_Spatial_DimPlot.png"),width=18,height=9, res = 400, units = "in")
p1 <- DimPlot(NPC, reduction = "umap", label = TRUE,pt.size = 1,
              # cols = cols,
              label.size = 10)
p2 <- SpatialDimPlot(NPC, label = TRUE, label.size = 8,cols = cols, pt.size.factor = 100)
print(p1 + p2)
dev.off()
