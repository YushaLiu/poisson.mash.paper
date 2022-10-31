# load in packages
devtools::load_all("poisson.mash.alpha")
setwd("applications_sc/Neutrophils")
library(Matrix)
library(scran)
library(glmpca)

### specifiy the selected conditions (i.e., cytokine treatments) measured on the same batch
trts <- c("Ctrl_2", "CCL20", "CXCL1", "CCL22", "CXCL5", "CCL11", "CCL4", "CCL17", "CCL5", "CXCL13", "CXCL10", "CXCL9",  
          "CXCL12", "GCSF", "MCSF", "GMCSF", "IFNg", "IL10", "IL12p70", "IL17a", "IL13", "IL15", "IL17f", "IL22",
          "IL18", "IL1a", "IL2", "IL3", "IL1b", "IL23", "IL21", "IL33", "IL25", "IL34", "IL36a", 
          "IL4", "IL6", "IL5", "IL7", "IL9", "IL11", "TGFb", "CCL2", "CCL3", "TSLP") 
R <- length(trts)

### load in the raw UMI count data and sample annotation information for all cell types
scdata <- readRDS("single-cell-cytokines/scdata.rds")
sample.info <- read.csv("single-cell-cytokines/whole_cyto_annot.csv")
sum(colnames(scdata)!=sample.info$X0)

### select the subset of Neutrophils in the selected treatment groups
idx.cell <- sample.info$cell_type=="Neutrophils" & sample.info$sample %in% trts
scdata <- scdata[, idx.cell]
sample.info <- sample.info[idx.cell,]

### remove genes with very low expression levels
idx.gene <- which(rowSums(scdata)>=25)
scdata <- scdata[idx.gene,]



################################ Aggregate UMI count data and size factors over conditions #########################################
### calculate the cell-wise library size using a deconvolution-based approach, which is more robust than taking the total UMI counts per cell
clusters <- quickCluster(scdata)
si <- calculateSumFactors(scdata, clusters=clusters)
names(si) <- colnames(scdata)
saveRDS(si, file = "size_factors.Rds")

### create a data object for poisson mash analysis, where cells from the same condition are aggregated to form "pseudo-bulk" data
conditions <- factor(sample.info$sample, levels=trts)
data <- pois_mash_set_data(scdata, conditions, si)
saveRDS(data, file = "data_jr.Rds")



############################## Perform GLM-PCA to estimate the matrix of latent factors for unwanted variation #################################
### construct the cell by condition design matrix which encodes the cell-wise assignment of treatment conditions
conditions <- factor(sample.info$sample, levels=trts)
data.X <- model.matrix(~conditions)
data.X <- data.X[,-1]
colnames(data.X) <- trts[-1]

### run glm-pca while adjusting for cell-specific size factors and gene-specific, condition-specific means
start_time = proc.time()
print("##########################################")
print("start fitting GLM-PCA for scdata from Neutrophils")
fit.glmpca <- glmpca(Y=scdata, X=data.X, L=4, fam="nb2", sz=si, ctl=list(verbose=TRUE, maxIter=600, tol=1e-6))
print("finish fitting GLM-PCA for scdata from Neutrophils")
runtime = proc.time() - start_time
fit.glmpca[["runtime"]] <- runtime
saveRDS(fit.glmpca, file = "glmpca.Rds")

### get the gene by factor matrix of latent factors causing unwanted variation, which is needed for subsequent poisson mash analysis
Fuv <- as.matrix(fit.glmpca$loadings)



################################ Run prefit step by ignoring fixed effects beta #########################################
start_time = proc.time()
print("##########################################")
print("start prefit without fixed effects for scdata from Neutrophils")
prefit <- pois_mash_ruv_prefit(data, Fuv, verbose=TRUE, control=list(maxiter=500))
print("finish prefit without fixed effects for scdata from Neutrophils")
runtime = proc.time() - start_time
prefit[["runtime"]] = runtime
saveRDS(prefit, file = "pois_mash_ruv_prefit.Rds")  



################################ Run ED step to estimate data-driven prior covariance matrices #########################################
### initialize data-driven prior covariance matrices
res.pca <- pois_cov_init(data, ruv=FALSE, npc=5, cutoff=abs(qnorm(0.05/2/R)))

### get the list of canonical prior covariances, which are by default rank-1 covariance matrices that model condition-specific effects
ulist.c <- pois_cov_canonical(data)

### combine all rank-1 prior covariances, and indicate which of them are data-driven ones
ulist <- c(res.pca$ulist, ulist.c)
ulist.dd <- c(rep(TRUE, length(res.pca$ulist)-1), rep(FALSE, R+1))

### run the ED step using only the subset of genes showing strong signals of differential expression
start_time = proc.time()
print("##########################################")
print("start fitting ED of possion mash ruv for scdata from Neutrophils")
fit.ed <- pois_cov_ed(data, subset=res.pca$subset, Ulist=res.pca$Ulist, ulist=ulist, ulist.dd=ulist.dd, ruv=TRUE, Fuv=Fuv, verbose=TRUE, 
                      control=list(maxiter=300, maxiter.q=25, maxpsi2=log(2), maxbias=1, tol.q=1e-2, tol.rho=1e-3, tol.stop=1e-6))
print("finish fitting ED of possion mash ruv for scdata from Neutrophils")
runtime = proc.time() - start_time
fit.ed[["runtime"]] = runtime
saveRDS(fit.ed, file = "pois_mash_ruv_ed.Rds")



################################ Run poisson mash RUV #########################################
### add a small epsilon to the diagonals of all data-driven prior covariance matrices, separately for matrices that are rank-1 or not
Ulist <- fit.ed$Ulist
H <- length(Ulist)
for(h in 1:H){
  Uh <- Ulist[[h]]
  Uh <- Uh/max(diag(Uh))
  Ulist[[h]] <- Uh + 1e-2*diag(R)
}
G <- length(fit.ed$ulist)
epsilon2.G <- rep(1e-8, G)
names(epsilon2.G) <- names(fit.ed$ulist)
epsilon2.G[ulist.dd] <- 1e-2

### run poisson mash RUV on all genes
start_time = proc.time()
print("##########################################")
print("start fitting poisson mash with ruv for scdata from Neutrophils")
res <- pois_mash(data=data, Ulist=Ulist, ulist=fit.ed$ulist, ulist.epsilon2=epsilon2.G, normalizeU=TRUE, gridmult=2, 
                 ruv=TRUE, Fuv=Fuv, rho=prefit$rho, update.rho = TRUE, verbose=TRUE, median_deviations = TRUE, 
                 init=list(mu=prefit$mu, psi2=prefit$psi2),
                 control=list(maxiter=300, maxiter.q=25, maxpsi2=log(2), maxbias=1, 
                              tol.q=1e-2, tol.rho=1e-3, tol.mu=1e-2, tol.psi2=2e-2, tol.bias=1e-2, nc=3)) 
print("finish fitting poisson mash with ruv for scdata from Neutrophils")
runtime = proc.time() - start_time
res[["runtime"]] = runtime
saveRDS(res, file = "pois_mash_ruv_fit_median.Rds")



print(sessionInfo())