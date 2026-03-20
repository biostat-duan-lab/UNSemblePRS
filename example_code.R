#####################################################
## Example code
## Implementation of UNSemblePRS
#####################################################

## Required R packages
library(kernlab)
library(sparsepca)


## Main function: UNSemblePRS()
## Input:
## X: n by p matrix of PRS predictions
## n: number of individuals
## p: number of PRS pre-trained models

UNSemblePRS <- function(X, standardize = TRUE, prop = 1){
  # Step 1: Scale each PRS to have zero mean and unit variance
  if(standardize){
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
  return(list(
    UNSemble = final.score, 
    weight = w.opt
    ))
}



#####################################################
## Simulation example
#####################################################

set.seed(123)

N <- 1000
p <- 100

## Generate genotype matrix
G_mat <- matrix(rbinom(N * p, 2, 0.3), N, p)

## True genetic effects
True_beta <- rnorm(p, 0, 1)

## True PRS
True_PRS <- G_mat %*% True_beta

## Simulate 5 PRS pre-trained models with different noise levels
PGS1_beta <- True_beta + rnorm(p, 0, 0.001)
PGS2_beta <- True_beta + rnorm(p, 0, 0.1)
PGS3_beta <- True_beta + rnorm(p, 0, 1)
PGS4_beta <- True_beta + rnorm(p, 0, 3)
PGS5_beta <- True_beta + rnorm(p, 0, 5)

## Compute predicted PRSs
PGS1 <- G_mat %*% PGS1_beta
PGS2 <- G_mat %*% PGS2_beta
PGS3 <- G_mat %*% PGS3_beta
PGS4 <- G_mat %*% PGS4_beta
PGS5 <- G_mat %*% PGS5_beta

## Apply UNSemblePRS
input <- cbind(PGS1, PGS2, PGS3, PGS4, PGS5)

output <- UNSemblePRS(input)

output$weight


## Evaluate performance
eval_mat <- cbind(
  True_PRS = True_PRS,
  UNSemble = output$UNSemble,
  PGS1 = PGS1,
  PGS2 = PGS2,
  PGS3 = PGS3,
  PGS4 = PGS4,
  PGS5 = PGS5
)

round(cor(eval_mat), 3)






