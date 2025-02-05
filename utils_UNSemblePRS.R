# Load required libraries
library(kernlab)
library(sparsepca)

# Input: p PRS model predictions for n individuals, denoted by matrix X (n x p)
UNSemblePRS <- function(X){
  # Step 1: Scale each PRS to have zero mean and unit variance
  X <- scale(X)
  
  # Step 2: Compute the feature space with the kernel mapping
  # Default kernel is the polynomial kernel of degree 3
  ## Polynomial kernel: polydot(degree = 3)
  ## Gaussian radial basis: rbfdot(sigma = 1)
  features <- kernelMatrix(polydot(degree = 3), x = t(X))
  
  # Step 3: Perform sparse PCA to obtain the aggregation weights
  # alpha.opt can be chosen by prior knowledge, 
  ## Default value: select half of the PRS models from PGS Catalog
  alpha.opt <- uniroot(f = function(alpha.opt) {
    w.opt <- rspca(X = features, alpha = alpha.opt, 
                   k = 1, verbose = FALSE)$loadings
    mean(w.opt != 0) - 0.5
  }, interval = c(0, 1))$root
  # Compute the optimized weights with the selected alpha.opt
  w.opt <- rspca(X = features, alpha = alpha.opt, 
                 k = 1, verbose = FALSE)$loadings
  
  # Step 4: Compute the final composite PRS score
  final.score <- X %*% w.opt
  
  # Output: The final aggregated PRS score
  return(final.score)
}
