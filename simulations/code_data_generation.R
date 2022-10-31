rm(list=ls())
setwd("simulations_sc")
library(Matrix)
library(seqgendiff)

################################ Read in data #########################################
scdata <- readRDS("raw_count.rds")
sample.info <- readRDS("sample_info.rds")

### specifiy the selected conditions
trts <- c("CCL11", "CCL17", "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CXCL5", "IFNg", "IL12p70", "IL18",
          "IL1a", "IL1b", "IL2", "IL3", "IL4", "IL5", "IL6", "IL7", "IL9", "IL27", "IL34", "IL36a", "TGFb", "TSLP") 

### take the subset of cells from given conditions
idx.cell <- sample.info$sample %in% trts
scdata <- scdata[, idx.cell]
sample.info <- sample.info[idx.cell,]
sum(colnames(scdata)!=sample.info$X0)
saveRDS(sample.info, file="data/sample_info.Rds")

### construct the design matrix without intercept
conditions <- factor(droplevels(sample.info$sample), levels=trts)
designmat <- model.matrix(~0+conditions)

### remove genes with totol counts over conditions fewer than 25 before thinning
keep.idx <- which(rowSums(scdata)>=25)
scdata <- scdata[keep.idx,]



################################ Specify treatment effects #########################################
### load in the prior covariances which are similar to those of the original data
ulist.sim <- readRDS("ulist_sim.rds")
ulist <- ulist.sim$ulist
ulist[[3]] <- pmax(ulist[[3]], 0)
pi.u <- rep(1/3, 3)

### determine distribution of effect sizes
wlist.upr <- log2(c(1.5, 2, 2.5, 3))
pi.w.upr <- c(0.2, 0.4, 0.3, 0.1)
wlist.lwr <- log2(c(2, 3, 4, 5))
pi.w.lwr <- rep(1/4, 4)

### pick the genes with relatively high or low expression to add fixed effects
gene.idx.upr <- which(rowSums(scdata) >= 400)
gene.idx.lwr <- which(rowSums(scdata) < 400 & rowSums(scdata) >= 200)



################################ Simulate single cell data #########################################
J <- nrow(scdata)
R <- length(trts)

### randomly permute the data and then add random effects
set.seed(2)

for (iter in 1:20){
  ### create null data by random permutation of cells with respect to the condition labels
  idx.perm <- sample(1:ncol(scdata), ncol(scdata), replace=FALSE)
  scdata.perm <- scdata[, idx.perm]
  
  ### simulate treatment effects
  beta <- matrix(0, nrow=J, ncol=R)
  
  non.null.idx.upr <- sort(sample(gene.idx.upr, 300, replace = FALSE))
  names(non.null.idx.upr) <- rownames(scdata)[non.null.idx.upr]
  non.null.idx.lwr <- sort(sample(gene.idx.lwr, 300, replace = FALSE))
  names(non.null.idx.lwr) <- rownames(scdata)[non.null.idx.lwr]
  non.null.idx <- c(non.null.idx.upr, non.null.idx.lwr)
  names(non.null.idx) <- rownames(scdata)[non.null.idx]
  
  for(j in 1:J){
    if(j %in% non.null.idx.upr){
      # simulate effect size
      w.j <- wlist.upr[which(as.numeric(rmultinom(1, 1, pi.w.upr))==1)]
      # add small noise to the effect size
      w.j <- rnorm(1, w.j, 0.1)
      
      # simulate effect-sharing pattern
      u.j <- ulist[[which(as.numeric(rmultinom(1, 1, pi.u))==1)]]
      beta[j,] <- w.j*ifelse(runif(1) > 0.5, 1, -1)*u.j
    }
    
    else if(j %in% non.null.idx.lwr){
      # simulate effect size
      w.j <- wlist.lwr[which(as.numeric(rmultinom(1, 1, pi.w.lwr))==1)]
      # add small noise to the effect size
      w.j <- rnorm(1, w.j, 0.1)
      
      # simulate effect-sharing pattern
      u.j <- ulist[[which(as.numeric(rmultinom(1, 1, pi.u))==1)]]
      beta[j,] <- w.j*ifelse(runif(1) > 0.5, 1, -1)*u.j      
    }
  }
  
  rownames(beta) <- rownames(scdata)
  colnames(beta) <- trts
  
  ### create signals through thinning
  coef_thin <- beta
  rownames(coef_thin) <- rownames(scdata)
  colnames(coef_thin) <- trts
  scdata.mod <- thin_diff(mat = as.matrix(scdata.perm), design_fixed = designmat, coef_fixed = coef_thin)
  X <- as(scdata.mod$mat, "sparseMatrix") 
  rownames(X) <- rownames(scdata)
  colnames(X) <- colnames(scdata)
  
  data <- list(X=X, beta=beta, non.null.idx=non.null.idx, non.null.idx.lwr=non.null.idx.lwr, non.null.idx.upr=non.null.idx.upr)

  saveRDS(data, file=paste0("data/raw_data", iter, ".Rds"))
}