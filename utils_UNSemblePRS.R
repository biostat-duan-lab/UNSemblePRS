# Main function
# Input: p PRS model predictions for n individuals, denoted by matrix X (n x p)
UNSemblePRS <- function(X, scale = TRUE, prop = 1){
  # Step 1: Scale each PRS to have zero mean and unit variance
  if(scale){
    X <- scale(X)
  }
  # Step 2: Compute the feature space with the kernel mapping
  # Default kernel is the polynomial kernel of degree 3
  ## Polynomial kernel: polydot(degree = 3)
  ## Gaussian radial basis: rbfdot(sigma = 1)
  features <- kernelMatrix(polydot(degree = 3), x = t(X))
  
  # Step 3: Perform sparse PCA to obtain the aggregation weights
  # alpha.opt can be chosen by prior knowledge, 
  ## Default value: select all of the PRS models
  alpha.opt <- uniroot(f = function(alpha.opt) {
    w.opt <- rspca(X = features, alpha = alpha.opt, 
                   k = 1, verbose = FALSE)$loadings
    mean(w.opt != 0) - prop
  }, interval = c(0, 1))$root
  # Compute the optimized weights with the selected alpha.opt
  w.opt <- rspca(X = features, alpha = alpha.opt, 
                 k = 1, verbose = FALSE)$loadings
  
  # Step 4: Compute the final composite PRS score
  final.score <- X %*% w.opt
  
  # Output: The final aggregated PRS score
  return(final.score)
}

# return the partial R2 of PRS scores given the age, sex and top 30 PCs
partial.R2 <- function(df, ypred,
                       continuous = TRUE,
                       cov.set = c('age', 'sex')){
  if(all(is.na(ypred))){
    # no information
    return(0)
  }else{
    cov.set <- paste0(cov.set, collapse = '+')
    # full model
    df_fit <- cbind(df, ypred = as.vector(ypred))
    if(continuous){
      # should place y pred to the end of linear regression
      lm.full.obj <- lm(
        paste0('pheno ~ ', cov.set, '+',
               paste0('PC', 1:30, collapse = '+'), '+ypred'),
        data = df_fit
      )
      # delete all NA
      df_noNA <- df %>% na.omit()
      # compute the sum of square in total
      SST <- sum((df_noNA$pheno - mean(df_noNA$pheno))**2)
      # compute the sequential increase in explained sum of squares
      anova(lm.full.obj)['ypred', ]$`Mean Sq`/SST 
      
    }else{
      # compute the adjusted AUC for binary outcome
      ypred.lm.obj <- lm(paste('ypred ~ ', cov.set, '+',
                               paste0('PC', 1:30, collapse = '+')),
                         data = df_fit)
      # obtain the residualized score
      ypred_residual <- as.vector(ypred) - predict(ypred.lm.obj,
                                                   newdata = df_fit)
      
      roc.binary(status = 'pheno',
                 variable = 'ypred_residual',
                 confounders = paste('~', cov.set, '+',
                                     paste0('PC', 1:30, collapse = '+')),
                 data = cbind(df_fit, 
                              ypred_residual = as.vector(ypred_residual)),
                 precision = seq(0.1, 0.9, by = 0.1))$auc
    }
  }
}

