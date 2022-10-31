setwd("simulations_sc")
library(Matrix)

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 
R <- length(trts)
sample.info <- readRDS("data/sample_info.Rds")
conditions <- factor(droplevels(sample.info$sample), levels=trts)


### apply Kruskal test 
for(idx in 1:20){
  data <- readRDS(file=paste0("data/raw_data", idx, ".Rds"))
  scdata <- data$X
  sum(colnames(scdata)!=sample.info$X0)
  
  start_time = proc.time()
  print("##########################################")
  print(sprintf("start running Kruskal test for replicate %d", idx))
  pval <- rep(NA, nrow(scdata))
  for(j in 1:nrow(scdata)){
    pval[j] <- kruskal.test(scdata[j,], conditions)$p.value
  }
  padj <- p.adjust(pval, method="fdr")
  print(sprintf("finish running Kruskal for replicate %d", idx))
  runtime = proc.time() - start_time
  saveRDS(list(df = data.frame(pval = pval, padj = padj, row.names = rownames(scdata)), runtime=runtime), 
          file = paste0("output/KM_test_rep_", idx, ".Rds"))
}



print(sessionInfo())