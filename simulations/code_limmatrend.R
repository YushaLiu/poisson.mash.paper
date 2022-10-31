setwd("simulations_sc")
library(Matrix)
library(edgeR)
library(limma)

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 
R <- length(trts)
sample.info <- readRDS("data/sample_info.Rds")
conditions <- factor(droplevels(sample.info$sample), levels=trts)
contrast.mat <- mashr::contrast_matrix(length(trts), ref = "CCL11", name = trts)  
contrast.mat <- t(contrast.mat)
colnames(contrast.mat) <- paste0(trts[-1], "-CCL11")


### apply limma trend
for(idx in 1:20){
  ### load in the data
  data <- readRDS(file=paste0("data/raw_data", idx, ".Rds"))
  scdata <- data$X
  sum(colnames(scdata)!=sample.info$X0)
  
  ### test if a gene is DE
  start_time = proc.time()
  dge <- DGEList(counts=scdata, group=conditions)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~conditions)
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  fit <- lmFit(y, design = design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- topTable(fit, n = Inf, adjust.method = "BH")
  runtime = proc.time() - start_time
  saveRDS(list(df = data.frame(pval = tt$P.Value, padj = tt$adj.P.Val, row.names = rownames(tt)), runtime=runtime), 
          file = paste0("output/limma_trend_rep_", idx, ".Rds"))
  
  ### test if a gene-condition is DE and estimate log fold change
  start_time = proc.time()
  dge <- DGEList(counts=scdata, group=conditions)
  dge <- calcNormFactors(dge)
  design <- model.matrix(~0+conditions)
  y <- new("EList")
  y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  fit <- lmFit(y, design = design)
  lmout.contrast <- limma::contrasts.fit(fit, contrast.mat)
  eout.contrast <- limma::eBayes(lmout.contrast, trend = TRUE, robust = TRUE)
  runtime = proc.time() - start_time
  saveRDS(list(bhat=eout.contrast$coef, pval=eout.contrast$p.value), file = paste0("output/limma_trend_ctrl_rep_", idx, ".Rds"))
}
