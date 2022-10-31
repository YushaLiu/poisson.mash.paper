args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

setwd("simulations_sc")
library(Matrix)
library(edgeR)

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 
R <- length(trts)
sample.info <- readRDS("data/sample_info.Rds")
conditions <- factor(droplevels(sample.info$sample), levels=trts)

### get the contrast matrix C to compare each treatment condition and the control
C <- mashr::contrast_matrix(length(trts), ref = "CCL11", name = trts) 

### load in data
data <- readRDS(file=paste0("data/raw_data", iter, ".Rds"))
scdata <- data$X
sum(colnames(scdata)!=sample.info$X0)

### test if a gene is DE
start_time = proc.time()
dge <- DGEList(counts=scdata, group=conditions)
dge <- calcNormFactors(dge)
cdr <- scale(colMeans(scdata > 0))
design <- model.matrix(~ cdr + conditions)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
qlf <- glmQLFTest(fit, coef = 3:ncol(design))
tt <- topTags(qlf, n = Inf)
runtime = proc.time() - start_time
saveRDS(list(tt = tt, df = data.frame(pval = tt$table$PValue, padj = tt$table$FDR, row.names = rownames(tt$table)), runtime=runtime), 
        file = paste0("output/edgeRQLFDetRate_rep_", iter, ".Rds"))

### test if a gene-condition is DE and estimate log fold change
start_time = proc.time()
dge <- DGEList(counts=scdata, group=conditions)
dge <- calcNormFactors(dge)
cdr <- scale(colMeans(scdata > 0))
design <- model.matrix(~ 0 + cdr + conditions)
dge <- estimateDisp(dge, design = design)
fit <- glmQLFit(dge, design = design)
contrast.mat <- rbind(0, t(C))
colnames(contrast.mat) <- rownames(C)

res <- list(NULL)
for(r in 1:ncol(contrast.mat)){
  qlf <- glmQLFTest(fit, contrast=contrast.mat[,r])
  res[[r]] <- qlf$table
}
names(res) <- colnames(contrast.mat)

runtime = proc.time() - start_time
saveRDS(list(res=res, runtime=runtime), file = paste0("output/edgeRQLFDetRate_ctrl_rep_", iter, ".Rds"))