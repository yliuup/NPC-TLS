library(monocle3)
library(Seurat)
library(dplyr)

load("cd8_filt.RData")

cd8_filt_sample <- subset(cd8_filt, cells=WhichCells(cd8_filt, downsample=2000))

subset_cells <- cd8_filt_sample
expression_matrix <- subset_cells@assays$RNA@counts %>% as.matrix()
cell_metadata <- subset_cells@meta.data
gene_annotation <- data.frame(gene_short_name = row.names(subset_cells@assays$RNA@counts), row.names = row.names(subset_cells@assays$RNA@counts), stringsAsFactors = F)


cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)


for (i in c(10,15,20,25,40,30,35,45)){
  # i = 25
  cds <- preprocess_cds(cds, num_dim = i)
  cds <- align_cds(cds, alignment_group = "Patient"
                   # ,residual_model_formula_str = "~ nCount_RNA + percent.mito"
  )
  
  
  cds <- reduce_dimension(cds,reduction_method = "UMAP",cores = 48,umap.fast_sgd = FALSE)
  
  
  
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  
  png(paste0("monocle3_dims_new",i,".png"),width=8,height=6,res = 400, units = "in")
  p <- plot_cells(cds,
                  color_cells_by = "annotation",
                  cell_size = 0.5,
                  group_label_size = 3,
                  label_groups_by_cluster=F,
                  label_leaves=F,
                  label_branch_points=F,
                  label_cell_groups=F)
  print(p)
  dev.off()
}

cols = c("#FED439FF","#FD7446FF","#C80813FF","#FD8CC1FF","#802268FF",
         "#D5E4A2FF", "#7E6148FF","#91331FFF","#F39B7FFF","#7F7F7FFF", "#D2AF81FF",
         "#370335FF","#FF7F0EFF","#A20056FF","#D6D6CEFF","#1B1919FF"
)

png(paste0("monocle3_dims_new",".png"),width=16,height=15,res = 400, units = "in")
plot_cells(cd8_cds,
           color_cells_by = "annotation",
           cell_size = 0.8,
           group_label_size = 3,
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F,
           label_cell_groups=F)+
  scale_color_manual(values = cols)

dev.off()


cds <- order_cells(cds)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time_bin="CD8_cluster"){
  cell_ids <- which(colData(cds)[, "annotation"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))

png("monocle3_pseudotime.png",width=8,height=6,res = 400, units = "in")
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=3,
           cell_size = 0.5)
dev.off()
