seurat_follow <- function(seurat_obj,
                          pathway = pathway,
                          dimss = dimss,
                          resolutions = resolutions ,
                          patientss = patientss,
                          regress_factor = regress_factor,
                          namess = namess){
  library(future)
  plan("multiprocess", workers = 10)
  plan()
  options(future.globals.maxSize = 10000000000000)
  
  TenXdat <- seurat_obj
  
  # TenXdat[["percent.mito"]] <- PercentageFeatureSet(TenXdat, pattern = "^MT-")
  
  TenXdat <- NormalizeData(object = TenXdat, normalization.method = "LogNormalize", scale.factor = 1e4)
  TenXdat <- FindVariableFeatures(object = TenXdat, selection.method = 'mean.var.plot', mean.cutoff = c(0.1, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 2000)
  TenXdat <- ScaleData(object = TenXdat, features = rownames(x = TenXdat), vars.to.regress = regress_factor)
  
  TenXdat <- RunPCA(object = TenXdat, features = VariableFeatures(object = TenXdat), verbose = FALSE)
  
  library(harmony)
  setwd(pathway)
  png("subset_cells_RunHarmony.png", width = 10, height = 10, res = 400, units = "in")
  TenXdat <- RunHarmony(TenXdat,patientss, plot_convergence = TRUE)
  dev.off()
  
  TenXdat <- FindNeighbors(TenXdat, reduction = "harmony", verbose = FALSE, dims = 1:dimss)
  
  TenXdat <- FindClusters(TenXdat, resolution = resolutions, verbose = FALSE, random.seed = 88)
  
  # plan("multiprocess", workers = 1) ####------------will be wrong, if workers != 1------------------####
  TenXdat <- RunUMAP(TenXdat, dims = 1:dimss, umap.method = "umap-learn",reduction = "harmony",
                     n.neighbors = 25L, min.dist = 0.5)
  
  names <- namess
  
  library(ggsci)
  png(paste0("dimplot_",names,"_ident.png"), width = 17, height = 15, res = 400, units = "in")
  p <- DimPlot(TenXdat,
          group.by = "ident", 
          pt.size = 1, 
          label = T,
          label.size = 10,
          raster = F,
          cols =  c(pal_npg()(10),pal_igv()(9),pal_uchicago("light")(9),pal_futurama()(12), pal_aaas()(10))[-8]
  )
  print(p)
  dev.off()
  return(TenXdat)
}
