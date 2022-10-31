args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

setwd("simulations_sc")
library(Matrix)
library(limma)
library(mashr)

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 
R <- length(trts)
sample.info <- readRDS("data/sample_info.Rds")
conditions <- factor(as.character(sample.info$sample), levels=trts)


################################ Read in data #########################################
data <- readRDS(file=paste0("data/raw_data", iter, ".Rds"))
scdata <- data$X
sum(colnames(scdata)!=sample.info$X0)

### log transform and normalize single cell count data
scale.factor <- colSums(scdata)
data.mash <- log(t(t(scdata)/scale.factor)*median(scale.factor) + 0.1)

# extract Bhat and Shat using limma  
data.X = model.matrix(~0+conditions)
colnames(data.X) <- trts
cov_of_interest = 1:ncol(data.X)
lmout = limma::lmFit(object = data.mash, design = data.X)
eout = limma::eBayes(lmout)
bhat = lmout$coefficients[,cov_of_interest,drop=FALSE]
shat = lmout$stdev.unscaled[,cov_of_interest,drop=FALSE] * sqrt(eout$s2.post)
rownames(shat) <- rownames(bhat)


################################ Run mash with common baseline ##################################
start_time = proc.time()
print("##########################################")
print(sprintf("start mash common baseline for replicate %d", iter))

### set up for mash 
mash.data = mash_set_data(bhat, shat, alpha = 0)
mash.data.L = mash_update_data(mash.data, ref = "CCL11")
mash.data.L = mash_set_data(Bhat=mash.data.L$Bhat, Shat=mash.data.L$Shat, alpha=1, V=mash.data.L$LSVSLt)

### canonical covariances
U.c = cov_canonical(mash.data.L)

### data-driven covariances
m.1by1 = mash_1by1(mash.data.L, alpha=1)
strong = get_significant_results(m.1by1, 0.05)
U.pca = cov_pca(mash.data.L, 5, subset=strong)
U.ed = cov_ed(mash.data.L, U.pca, subset=strong) 
for(k in 1:length(U.ed)){
  Uk <- U.ed[[k]]
  Uk <- max(diag(Uk))*(Uk/max(diag(Uk)) + 1e-2*diag(ncol(Uk)))
  U.ed[[k]] <- (Uk + t(Uk))/2
}

### fit mash
mash.fit = mash(mash.data.L, c(U.c,U.ed), algorithm.version = 'R')
print(sprintf("finish mash common baseline for replicate %d", iter))
runtime = proc.time() - start_time
mash.fit[["runtime"]] = runtime
saveRDS(mash.fit, file = paste0("output/mash_cb_ctrl_rep_", iter, ".Rds"))  



print(sessionInfo())