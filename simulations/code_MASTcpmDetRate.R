args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

setwd("simulations_sc")
library(Matrix)
library(edgeR)
library(MAST)

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 
R <- length(trts)
sample.info <- readRDS("data/sample_info.Rds")
conditions <- factor(droplevels(sample.info$sample), levels=trts)

### load in the data
data <- readRDS(file=paste0("data/raw_data", iter, ".Rds"))
scdata <- data$X
sum(colnames(scdata)!=sample.info$X0)

### test if a gene is DE
start_time = proc.time()
cdr <- scale(colMeans(scdata > 0))
dge <- DGEList(counts = scdata)
dge <- edgeR::calcNormFactors(dge)
sca <- FromMatrix(exprsArray = edgeR::cpm(dge, log=TRUE, prior.count=1), cData = data.frame(wellKey = colnames(scdata), grp = conditions, cdr = cdr))
zlmCond <- zlm(~ cdr + grp, sca) 
mast <- lrTest(zlmCond, "grp")
runtime = proc.time() - start_time
saveRDS(list(res = mast, df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"], row.names = names(mast[, "hurdle", "Pr(>Chisq)"])), runtime=runtime), 
        file = paste0("output/MASTcpmDetRate_rep_", iter, ".Rds"))

### test if a gene-condition is DE and estimate log fold change
start_time = proc.time()
cdr <- scale(colMeans(scdata > 0))
dge <- DGEList(counts = scdata)
dge <- edgeR::calcNormFactors(dge)
sca <- FromMatrix(exprsArray = edgeR::cpm(dge, log=TRUE, prior.count=1), 
                  cData = data.frame(wellKey = colnames(scdata), grp = conditions, cdr = cdr))
zlmCond <- zlm(~ 0 + cdr + grp, sca) 
contrast1 <- rbind(0, diag(R))
contrast1 <- contrast1[, -c(1)]
colnames(contrast1) <- paste0(trts[-1], "-", trts[1])
lfc <- getLogFC(zlmCond, contrast0=c(0, 1, rep(0, R-1)), contrast1=contrast1)
runtime = proc.time() - start_time
saveRDS(list(lfc = lfc, runtime = runtime), file = paste0("output/MASTcpmDetRate_ctrl_rep_", iter, ".Rds"))



print(sessionInfo())