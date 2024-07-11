library(dplyr)
library(Seurat)
library(CellChat)


###Pre-process
load("CellChat_CD4_B_t_raw.RData")

cellchat_t <- createCellChat(object = data.input, meta = identity, group.by = "group")
##
cellchat_t <- addMeta(cellchat_t, meta = identity, meta.name = "labels")
cellchat_t <- setIdent(cellchat_t, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_t@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat_t@idents)) # number of cells in each cell group
##


# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
# cellchat@DB <- CellChatDB.use

cellchat_t@DB <- CellChatDB.human


load("CellChat_CD4_B_n_raw.RData")
cellchat_n <- createCellChat(object = data.input, meta = identity, group.by = "group")

cellchat_n <- addMeta(cellchat_n, meta = identity, meta.name = "labels")
cellchat_n <- setIdent(cellchat_n, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_n@idents) # show factor levels of the cell labels

groupSize <- as.numeric(table(cellchat_n@idents)) # number of cells in each cell group
##


# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
# cellchat@DB <- CellChatDB.use

cellchat_n@DB <- CellChatDB.human

###work_follow
cellchat_follow <- function(cellchat = cellchat){
  cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
  library(future)
  future::plan("multiprocess", workers = 10) # do parallel  
  options(future.globals.maxSize = 10000000000000)
  
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)  
  
  ###
  cellchat <- computeCommunProb(cellchat,population.size = TRUE,type = "truncatedMean" , trim = 0.1) 
  # mycomputeCommunProb <-edit(computeCommunProb)  
  # environment(mycomputeCommunProb) <- environment(computeCommunProb)
  # cellchat <- mycomputeCommunProb(cellchat)  
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
  
  
}
cellchat_t <- cellchat_follow(cellchat_t)
cellchat_n <- cellchat_follow(cellchat_n)

save(cellchat_t,file = "cellchat_t.RData")
save(cellchat_n,file = "cellchat_n.RData")




###Merge
load("cellchat_t.RData")
load("cellchat_n.RData")
object.list <- list(Tumour = cellchat_t, Normal = cellchat_n)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


###Plot
pdf('CellChat_compareInteractions.pdf',width = 20,height =10)
#
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
dev.off()


pdf('CellChat_rankNet.pdf',width = 20,height =15)

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()


pdf('CellChat_netVisual_bubble.pdf',width = 10,height =6)

gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 2, 
                        title.name = "Increased signaling in LS", angle.x = 45, remove.isolate = T)

gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 1, 
                        title.name = "Decreased signaling in LS", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2

dev.off()



# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Tumour"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.01, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "Tumour",ligand.logFC = 0.01, receptor.logFC = 0.01)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "Normal",ligand.logFC = -0.01, receptor.logFC = -0.01)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pdf('CellChat_netVisual_bubble_DEG_normal.pdf',width = 20,height =15)


pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.up$interaction_name <- as.character(pairLR.use.up$interaction_name)
# pairLR.use.up_filt <- pairLR.use.up[grep(pairLR.use.up$interaction_name,pattern="ITG",invert = T)]

pairLR.use.up_filt <- pairLR.use.up[grep(pairLR.use.up$interaction_name,pattern="ITG|COL|EFNB|LAM|MIF",invert = T),]  %>% as.data.frame
colnames(pairLR.use.up_filt) <- "interaction_name"


gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1,2), comparison = c(1, 2), 
                        angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]

pairLR.use.down$interaction_name <- as.character(pairLR.use.down$interaction_name)

pairLR.use.down_filt <- pairLR.use.down[grep(pairLR.use.down$interaction_name,pattern="ITG|COL|EFNB|LAM|MIF",invert = T),]  %>% as.data.frame
colnames(pairLR.use.down_filt) <- "interaction_name"


gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(3), targets.use = c(1,2), comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()
