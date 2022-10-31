args = commandArgs(trailingOnly=TRUE)
iter = as.integer(args[1])

devtools::load_all("poisson.mash.alpha")
setwd("simulations_sc")
library(Matrix)
library(scran)    
library(glmpca)


################################ Read in functions and data #########################################
### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 
R <- length(trts)
sample.info <- readRDS("data/sample_info.Rds")

### get the contrast matrix C to compare each treatment condition and the control
C <- mashr::contrast_matrix(length(trts), ref = "CCL11", name = trts)   

# construct the vector of conditions
conditions <- factor(droplevels(sample.info$sample), levels=trts)

# determine the value of epsilon to add to the diagonal of prior covariances
epsilon <- 1e-2

### read in a simulated single cell dataset
data <- readRDS(file=paste0("data/raw_data", iter, ".Rds"))
scdata <- list(Y=data$X, condition=conditions)

### compute cell-specific size factors using scrna
clusters <- quickCluster(scdata$Y)
si <- calculateSumFactors(scdata$Y, clusters=clusters)

### create a data object for poisson mash analysis
data <- pois_mash_set_data(scdata$Y, conditions, si)



############################ Estimate the matrix of latent factors causing unwanted variation #####################################
### run glmpca to estimate matrix of latent factors while adjusting for cell-specific size factors and gene-specific, condition-specific means
print("##########################################")
print("start fitting GLMPCA to estimate matrix of latent factors")
design <- model.matrix(~scdata$condition)
fit.glmpca <- glmpca(Y=scdata$Y, X=design[,-1], L=4, fam="nb2", sz=si, ctl=list(verbose=TRUE, tol=1e-5))
Fuv <- as.matrix(fit.glmpca$loadings)
print("finish fitting GLMPCA to estimate matrix of latent factors")
saveRDS(fit.glmpca, file = paste0("output/glmpca_rep", iter, ".Rds"))



############################ Prefit the model to initialize parameters #####################################
print("##########################################")
print("start prefit for parameter initialization")
prefit <- pois_mash_ruv_prefit(data, Fuv, verbose=TRUE)
print("finish prefit for parameter initialization")



############################ Estimate data-driven prior covariance matrices #####################################
### initialize data-driven prior covariance matrices
res.pca <- pois_cov_init(data, ruv=TRUE, Fuv=Fuv, rho=prefit$rho, npc=5)

### combine all the rank-1 prior covariance matrices 
ulist.c <- pois_cov_canonical(data)
ulist <- c(res.pca$ulist, ulist.c)
ulist.dd <- c(rep(TRUE, length(res.pca$ulist)-1), rep(FALSE, R+1))

### run the ED step
print("##########################################")
print("start fitting ED step")
fit.ed <- pois_cov_ed(data, subset=res.pca$subset, Ulist=res.pca$Ulist, ulist=ulist, ulist.dd=ulist.dd,
                      ruv=TRUE, Fuv=Fuv, verbose=TRUE, init=prefit)
print("finish fitting ED step")



############################ Run Poisson mash RUV #####################################
# add epsilon*I to each full rank prior covariance matrix
Ulist <- fit.ed$Ulist
H <- length(Ulist)
for(h in 1:H){
  Uh <- Ulist[[h]]
  Uh <- Uh/max(diag(Uh))
  Ulist[[h]] <- Uh + epsilon*diag(R)
}

# add epsilon*I to rank-1 prior covariance matrix
G <- length(fit.ed$ulist)
epsilon2.G <- rep(1e-8, G)
names(epsilon2.G) <- names(fit.ed$ulist)
epsilon2.G[ulist.dd] <- epsilon

start_time = proc.time()
print("##########################################")
print("start fitting poisson mash with ruv")
res <- pois_mash(data=data, Ulist=Ulist, ulist=fit.ed$ulist, ulist.epsilon2=epsilon2.G, normalizeU=TRUE, gridmult=2.5, 
                 ruv=TRUE, Fuv=Fuv, rho=prefit$rho, verbose=TRUE, C=C, res.colnames=rownames(C),
                 init=list(mu=prefit$mu, psi2=prefit$psi2), control=list(maxiter=300, tol.mu=1e-2, tol.psi2=2e-2, tol.bias=1e-2))
print("finish fitting poisson mash with ruv")
runtime = proc.time() - start_time
res[["runtime"]] = runtime
saveRDS(res, file = paste0("output/pois_mash_ruv_fit_ctrl_rep", iter, ".Rds"))



print(sessionInfo())