setwd("applications_sc/Neutrophils")
library(Matrix)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(ggrepel)

### load in the color code for conditions
cols.trt <- readRDS("condition_colors.Rds")


############################### Look at the poisson mash ruv fit to neutrophils data ################################
### load in the fit from ED
fit.ed <- readRDS("pois_mash_ruv_ed.Rds")

### check ELBO for ED
plot(fit.ed$ELBO, xlab = "niter", ylab = "elbo")

### load in the entire fit
res <- readRDS("pois_mash_ruv_fit_median.Rds")
post <- res$result$beta_median_dev_post
sum(names(cols.trt)!=colnames(post$PosteriorMean))

### check ELBO for the entire fit
plot(res$pois.mash.fit$ELBO[res$pois.mash.fit$ELBO!=0], xlab = "niter", ylab = "elbo")

### check estimated weights of the mash prior
pi.mat <- res$pois.mash.fit$pi
pheatmap(pi.mat, cluster_rows=FALSE, cluster_cols=FALSE, fontsize_row = 6, fontsize_col=6, 
         main="Estimated weights of prior covariances in poisson mash ruv (Neutrophils)")
rowSums(pi.mat)


### plot the top eigenvectors of the single data-driven prior covariance matrix with the largest estimated weight
eig.tPCA <- eigen(fit.ed$Ulist[[1]])
pve.tPCA <- eig.tPCA$values/sum(eig.tPCA$values)

pdf("prior_covariances.pdf", width=12, height=12)
par(mfrow=c(2,2))
for (k in 1:4){
  v <- eig.tPCA$vectors[,k]
  barplot(v/v[which.max(abs(v))], names = names(cols.trt), cex.names = 0.4, las = 2, col = cols.trt,
          main = paste0(names(res$pois.mash.fit$Ulist)[1], ": weight ", round(rowSums(pi.mat)[1],3),  
                        ":\nEigenvector ", k, " (pve = ", round(pve.tPCA[k],3), ")"))
}
dev.off()


### look at the posterior mean and standard deviations of DE effects (relative to the median) for all genes
pm <- post$PosteriorMean
psd <- post$PosteriorSD


### get the list of genes which are identified as DE in at least one conditions based on local false sign rate
lfsr <- ashr::compute_lfsr(post$NegativeProb, post$ZeroProb)
idx.pois <- which(apply(lfsr, 1, min) < 0.05)
length(idx.pois)
idx.pois <- idx.pois[order(apply(lfsr[idx.pois,], 1, min))]     # order identified DE genes by significance