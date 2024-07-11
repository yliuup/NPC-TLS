library(pROC)
library(survival)
library(survminer)
library(data.table)
library(ggplot2)
library(ggsignif)
library(ggpubr)

NPC_RNAseq <- read.table("TPM.txt", header = T, sep = "\t", stringsAsFactors = F, row.names = 1)


clinical_features <- read.table("SRR_survival_statue.txt", row.names = 1, header = T, sep = "\t", stringsAsFactors = F)

names(NPC_RNAseq) <- do.call(rbind, strsplit(names(NPC_RNAseq), split = "_"))[,1]
samples_of_RNA_seq <- names(NPC_RNAseq)

overlap_samples <- intersect(row.names(clinical_features), samples_of_RNA_seq)
Rna_seq <- NPC_RNAseq[, overlap_samples]

clinical_features_overlaps <- clinical_features[overlap_samples, ]
clinical_features_overlaps$XIST <- ""
clinical_features_overlaps$XIST <- as.numeric(Rna_seq["XIST", ])> 1

clinical_features_overlaps[clinical_features_overlaps$XIST, "XIST"] <- "Female"
clinical_features_overlaps[clinical_features_overlaps$XIST==FALSE, "XIST"] <- "Male"

strs <- c("Male", "Female")
names(strs) <- c("TRUE", "FALSE")


expression_dataset_zscore <- t(apply(Rna_seq[gene_used, ], 1, ZScore)) ## z-transform each genes
expression_dataset_zscore_mean <- apply(expression_dataset_zscore, 2, mean) ## average the expression of geneset

normalize.by.gene <- "PTPRC"
expression_dataset_zscore_mean <- expression_dataset_zscore_mean - as.vector(apply(Rna_seq[normalize.by.gene,], 1, ZScore))



roc.module <- roc(clinical_features_overlaps$DFS, as.numeric(expression_dataset_zscore_mean))

coords.modult <- coords(roc.module, "best", 
                        input=c("threshold"),
                        ret=c("threshold", "specificity", "sensitivity"),
                        as.list=FALSE, drop=TRUE, best.method=c("closest.topleft"))



binary_expression <- as.numeric((expression_dataset_zscore_mean > coords.modult["threshold"])+0)



clinical_features_overlaps$genes_exp <- binary_expression

km.as.gene <- with(clinical_features_overlaps, survfit(Surv(DFS_time, DFS) ~ genes_exp, data = clinical_features_overlaps, conf.type = "log-log"))

ggsurvplot(km.as.gene, conf.int=F, pval=TRUE, risk.table=TRUE,
           legend.labs=c("0", "1"), legend.title=genes,
           palette=c("blue", "red"),
           main="Kaplan-Meier Curve for NPC Survival",
           risk.table.height=.20)
