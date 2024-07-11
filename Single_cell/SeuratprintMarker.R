SeuratprintMarker <- function(pbmc10k,cores = 20 , fc.thresholds = fc.thresholds){
  plan("multiprocess", workers = 1)
  require(parallel)
  date()
  system.time({
    date()
    result <- mclapply(as.character(levels(pbmc10k@active.ident)),
                       FUN =  function(x) {FindMarkers(pbmc10k, ident.1 = x, ident.2 = NULL)},
                       mc.cores = cores)
    date()
  })
  date()
  RESULT <- result
  
  roundN <- 1
  while(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
    if(any(mapply(length, result, SIMPLIFY = TRUE)!=5)){
      recalculate_clusters <- which(mapply(length, result, SIMPLIFY = TRUE)!=5)-1
      print(recalculate_clusters)
      result1 <- mclapply(recalculate_clusters,
                          FUN =  function(x) {FindMarkers(pbmc10k, ident.1 = x, ident.2 = NULL)},
                          mc.cores = cores)
    }
    print(roundN + 1)
    for(i in 1:length(recalculate_clusters)){
      result[c(recalculate_clusters+1)[i]] <- result1[i]
    }
  }
  
  all_markers <- do.call(rbind, result)
  all_markers$gene <- unlist(mapply(rownames, result))
  all_markers$cluster <- rep(levels(pbmc10k@active.ident), times = mapply(dim, result, SIMPLIFY = TRUE)[1,])
  # return(all_markers)
  TOP_N <- function(x, n, pct.1 = 0.1,pct.2 = 1, first = "avg_log2FC", second = "p_val_adj", sig.padj = NULL, fc.threshold = c(-0.25, 0.25)){
    if(table(names(x) %in% c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "gene", "cluster"))["TRUE"] == 7){
      reordered_daf <- data.frame()
      if (!is.null(pct.1) & (pct.1 > 0 & pct.1 < 1)) {
        x <- x[x$pct.1 >= pct.1, ]
      }
      if (!is.null(pct.2) & (pct.2 > 0 & pct.2 < 1)) {
        x <- x[x$pct.2 <= pct.2, ]
      }
      if (!is.null(sig.padj) & is.numeric(sig.padj)) {
        x <- x[x$p_val_adj <= sig.padj, ]
      }
      if (length(fc.threshold) == 2) { 
        fc.threshold <- sort(fc.threshold)
        x <- x[(x$avg_log2FC < fc.threshold[1] | x$avg_log2FC > fc.threshold[2]), ]
      } else if (length(fc.threshold) == 1) {
        x <- x[x$avg_log2FC >= fc.threshold, ]
      }
      
      for (i in unique(x$cluster)) {
        if(n > dim(x[x$cluster == i, ])[1]){
          message("FBI warning, n < ", n)
          tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first), get(second)))[1:(dim(x[x$cluster == i, ])[1]), ]
          reordered_daf <- rbind(reordered_daf, tmp)
        } else {
          tmp <- dplyr::arrange(.data = x[x$cluster == i, ], desc(get(first), get(second)))[1:n, ]
          reordered_daf <- rbind(reordered_daf, tmp)
        }
      }
    }
    else {
      message("be careful !")
    }
    return(reordered_daf)
  }
  all_markers_sort <- all_markers %>% TOP_N(10000, pct.1 = 0.1, sig.padj= 0.05, fc.threshold = fc.thresholds)
  return(all_markers_sort)
}
