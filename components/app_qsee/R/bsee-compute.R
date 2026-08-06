
bsee_compute_batchcorrect <- function(X, samples, pheno, progress=NULL) {

  if(length(pheno)==1 && pheno %in% colnames(samples)) {
    pheno <- samples[,pheno]
  }
  pheno <- as.character(pheno)
  
  ## default normalization
  X <- playbase::counts.mergeDuplicateFeatures(X)
  X <- limma::normalizeQuantiles(playbase::logCPM(X))
  X <- playbase::svdImpute2(X) ## standard impute
  
  ## svd results
  res <- irlba::irlba(X, nv=4)
  colnames(res$v) <- paste0("PC",1:ncol(res$v))
  colnames(res$u) <- paste0("PC",1:ncol(res$u))
  rownames(res$u) <- rownames(X)
  rownames(res$v) <- colnames(X)
  
  if(!is.null(progress)) progress$set(message = paste("Analyzing PCA..."), value = 0.3)
  
  ## PC-traits correlation
  H <- playbase::expandPhenoMatrix(samples)
  R <- cor( H, res$v )
  P <- matrix(NA, ncol(H), ncol(res$v))
  colnames(P) <- colnames(res$v)
  rownames(P) <- colnames(H)
  i=1
  for(i in 1:ncol(H)) {
    m <- playbase::gx.limma( t(res$v), H[,i], fdr=1, lfc=0,
      sort.by='none')
    P[i,] <- m$P.Value
  }
  res$P <- P  
  res$X <- X
  res$samples <- samples
  
  analyze_batch_effects <- function(res) {    
    samples <- res$samples
    X <- res$X
    bc <- playbase::detectBatchEffects(
      X = X,
      samples = samples,
      pheno = pheno,
      ## contrasts = NULL,
      ## params = c("statistical","technical","pca"),
      params = c("statistical", "technical"),
      p.pca = 0.5,
      p.pheno = 0.05,
      k.pca = 10,
      nv = 2,
      xrank = NULL
    )
    return(bc)
  }
  
  compare_methods <- function(res) {    
    X0 <- res$X
    ntop_features <- 1000    
    methods <- c("ComBat", "limma", "RUV", "SVA", "NPM")
    ##if (ncol(X0) > 100) methods <- methods[methods != "NPM"]
    xlist.init <- list("uncorrected" = X0)   
    batch.pars <- unlist(res$effects$params)    
    mres <- playbase::compare_batchcorrection_methods(
      X0,
      res$samples,
      pheno = pheno,
      contrasts = NULL,
      batch.pars = batch.pars,
      clust.method = "pca",
      methods = methods,
      #evaluate = FALSE, ## no score computation
      evaluate = TRUE, ## no score computation
      xlist.init = xlist.init,
      ntop = ntop_features,
      npc = 5
    )    
    return(mres)
  }

  if(!is.null(progress))
    progress$set(message = paste("Analyzing batch effects..."), value = 0.5)
  res$effects <- analyze_batch_effects(res)

  if(!is.null(progress)) 
    progress$set(message = paste("Comparing methods..."), value = 0.8)
  res$methods <- compare_methods(res)

  message("[bsee_compute_batchcorrect] done")
  return(res)
}
