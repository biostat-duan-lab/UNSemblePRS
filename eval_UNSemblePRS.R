# evaluating the performance of PRS model for each trait
# load required libraries
library(RISCA)
library(kernlab)
library(sparsepca)
library(doParallel)
library(irlba)
source('utils_UNSemblePRS.R')
# setting
niters <- 100 # number of replications
Nsub <- 1000 # size of random samples
kernel_fun <- polydot(degree = 3) # kernel function for feature mapping
# phenotype of interest
# continuous traits: height, BMI, HDL, LDL, Diastolic Blood Pressure (DBP)
# binary traits: breast cancer (BC), Inflammatory Bowel Disease (IBD),
## Chronic kidney disease (CKD), Atrial fibrillation (AF), 
## Coronary artery disease (CAD)
pheno <- 'height' 
if(pheno %in% c('height', 'bmi', 'hdl', 'ldl', 'dbp')){
  continuous <- TRUE
}
if(pheno %in% c('BC', 'IBD', 'CKD', 'AF', 'CAD')){
  continuous <- FALSE
}
# prepare the full PRS dataframe
## 1. person_id and the PRS scores from the PGS Catalog
## 2. phenotype of the trait of interest
## 3. Patient information: ancestry, age, gender and top genetic PCs
### age should be larger than 21,
### gender is the sex at birth
### an example of PRS dataset:
head(dfPRS)
#    person_id  PGS000297   PGS001229 pheno ancestry      age     sex             PC1
# 1   1000004 0.00833752 1.32662e-04 1.872      eur 83.54278      Male   -0.000174243
# 2   1000033 0.00814360 1.04787e-04 1.874      eur 67.54278      Male   -0.000171157
# 3   1000059 0.00848452 1.23336e-04 1.645      eur 66.54346    Female   -0.000169760
# 4   1000061 0.00801136 1.42889e-06 1.803      eur 59.54278      Male   -0.000172519
# 5   1000070 0.00829886 8.75520e-05 1.588      eur 76.54483 PMI: Skip   -0.000174638
# 6   1000079 0.00797628 1.10177e-04 1.884      eur 31.54278      Male   -0.000180965

n_cores <- detectCores()
# register the cluster
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)
# begin the evaluation in parallel 
r2.est <- foreach(seed = 1:niters,
                  .combine = 'rbind',
                  .packages = c("dplyr", "kernlab", 'RISCA', 'irlba', 'sparsepca')) %dopar% {
                    # select a subset from the full PRS dataframe
                    dfPRS_sample <- dfPRS[sample(1:nrow(dfPRS), 
                                                 size = Nsub, 
                                                 replace = TRUE), ]
                    # obtain the scaled PRS scores
                    PRS_X <- scale(dfPRS_sample %>% select(starts_with('PGS')))
                    
                    # obtain the best PRS model
                    ## Ancestry
                    bestAncestry <- with(dfPRS_sample,
                                         by(dfPRS_sample, ancestry, function(df){
                                           # compute the sequential R2 for each PGS model
                                           r2.list <- sapply(grep('PGS', colnames(df), 
                                                                  value = TRUE), function(pgs.name){
                                                                    partial.R2(df = df,
                                                                               ypred = df[, pgs.name],
                                                                               continuous = continuous)
                                                                  })
                                           # return the name of PRS model with the high partial R2
                                           (grep('PGS', colnames(df), value = TRUE))[r2.list %>% which.max()]
                                         })) %>% array2DF()
                    ## Age stratum
                    dfPRS_sample <- dfPRS_sample %>% mutate(AgeStrata = cut(age, 
                                                                            breaks = c(20, 40, 60, 80, Inf),
                                                                            labels = paste0('age', 1:4)))
                    bestAge <- with(dfPRS_sample,
                                    by(dfPRS_sample, AgeStrata, function(df){
                                      # compute the sequential R2 for each PGS model
                                      r2.list <- sapply(grep('PGS', colnames(df), 
                                                             value = TRUE), function(pgs.name){
                                                               partial.R2(df = df,
                                                                          ypred = df[, pgs.name],
                                                                          continuous = continuous)
                                                             })
                                      # return the name of PRS model with the high partial R2
                                      (grep('PGS', colnames(df), value = TRUE))[r2.list %>% which.max()]
                                    })) %>% array2DF()
                    ## sex
                    bestSex <- with(dfPRS_sample,
                                       by(dfPRS_sample, sex, function(df){
                                         # compute the sequential R2 for each PGS model
                                         r2.list <- sapply(grep('PGS', colnames(df), 
                                                                value = TRUE), function(pgs.name){
                                                                  partial.R2(df = df,
                                                                             ypred = df[, pgs.name],
                                                                             continuous = continuous,
                                                                             cov.set = c('age'))
                                                                })
                                         # return the name of PRS model with the high partial R2
                                         (grep('PGS', colnames(df), value = TRUE))[r2.list %>% which.max()]
                                       })) %>% array2DF()
                    
                    # evaluate the performance for each ancestry stratum
                    r2.estAncestry <- lapply(c('eur', 'afr', 'amr', 'eas'), function(ancestry){
                      
                      # 1) average of the PRS scores
                      r2.avg <- partial.R2(df = dfPRS_sample[dfPRS_sample$ancestry == ancestry,],
                                             ypred = PRS_X[dfPRS_sample$ancestry == ancestry, ]%>%
                                               apply(1, mean),
                                             continuous = continuous)
                      
                      # 2) aggregation of PRS scores by PCA weights
                      w <- abs(irlba(PRS_X[dfPRS_sample$ancestry == ancestry, ],
                                     nv = 1)$v)
                      r2.PCA <- partial.R2(df = dfPRS_sample[dfPRS_sample$ancestry == ancestry,],
                                             ypred = (PRS_X %*% w)[dfPRS_sample$ancestry == ancestry],
                                             continuous = continuous)
                      # 3) the best single PRS score
                      r2.best <- partial.R2(df = dfPRS_sample[dfPRS_sample$ancestry == ancestry,],
                                              ypred = PRS_X[dfPRS_sample$ancestry == ancestry,
                                                            bestAncestry$Value[which(bestAncestry$ancestry == ancestry)]],
                                              continuous = continuous)
                      # 4) the unsupervised ensemblePRS
                      s.UNSemblePRS <- UNSemblePRS(PRS_X[dfPRS_sample$ancestry == ancestry, ])
                      r2.UNSemblePRS <- partial.R2(df = dfPRS_sample[dfPRS_sample$ancestry == ancestry,],
                                                     ypred = s.UNSemblePRS,
                                                     continuous = continuous)
                      data.frame(context = ancestry,
                                 best = r2.best,
                                 avg = r2.avg,
                                 PCA = r2.PCA,
                                 UNsemblePRS = r2.UNSemblePRS)
                    }) %>% 
                      do.call(rbind, .)
                    
                    # evaluate the performance for each age stratum
                    r2.estAge <- lapply(c('age1', 'age2', 'age3', 'age4'), FUN = function(AgeStrata){
                      
                      if(AgeStrata == 'age1'){age_range <- c(20, 40)}
                      if(AgeStrata == 'age2'){age_range <- c(40, 60)}
                      if(AgeStrata == 'age3'){age_range <- c(60, 80)}
                      if(AgeStrata == 'age4'){age_range <- c(80, Inf)}
                      
                      # 1) average of the PRS scores
                      r2.avg <- partial.R2(df = dfPRS_sample[dfPRS_sample$age >= age_range[1] &
                                                                 dfPRS_sample$age < age_range[2], ],
                                             ypred = PRS_X[dfPRS_sample$age >= age_range[1] &
                                                             dfPRS_sample$age < age_range[2], ]%>%
                                               apply(1, mean),
                                             continuous = continuous)
                      # 2) aggregation of PRS scores by PCA weights
                      w <- abs(irlba(PRS_X[dfPRS_sample$age >= age_range[1] & 
                                             dfPRS_sample$age < age_range[2], ],
                                     nv = 1)$v)
                      r2.PCA <- partial.R2(df = dfPRS_sample[dfPRS_sample$age >= age_range[1] &
                                                                 dfPRS_sample$age < age_range[2], ],
                                             ypred = (PRS_X %*% w)[dfPRS_sample$age >= age_range[1] & 
                                                                     dfPRS_sample$age < age_range[2]],
                                             continuous = continuous)
                      # 3) the best single PRS score
                      r2.best <- partial.R2(df = dfPRS_sample[dfPRS_sample$age >= age_range[1] &
                                                                  dfPRS_sample$age < age_range[2], ],
                                              ypred = PRS_X[dfPRS_sample$age >= age_range[1] & 
                                                              dfPRS_sample$age < age_range[2],
                                                            bestAge$Value[which(bestAge$AgeStrata == AgeStrata)]],
                                              continuous = continuous)
                      # 4) the unsupervised ensemblePRS
                      s.UNSemblePRS <- UNSemblePRS(PRS_X[dfPRS_sample$age >= age_range[1] &
                                                           dfPRS_sample$age < age_range[2], ])
                      r2.UNSemblePRS <- partial.R2(df = dfPRS_sample[dfPRS_sample$age >= age_range[1] &
                                                                         dfPRS_sample$age < age_range[2], ],
                                                     ypred = s.UNSemblePRS,
                                                     continuous = continuous)
                      data.frame(context = AgeStrata,
                                 best = r2.best,
                                 avg = r2.avg,
                                 PCA = r2.PCA,
                                 UNsemblePRS = r2.UNSemblePRS)
                    }) %>% 
                      do.call(rbind, .)
                    
                    # evaluate the performance for each gender
                    r2.estSex <- lapply(c('Male', 'Female'), FUN = function(SexStrata){
                      
                      # 1) average of the PRS scores
                      r2.avg <- partial.R2(df = dfPRS_sample[dfPRS_sample$sex == SexStrata,],
                                             ypred = PRS_X[dfPRS_sample$sex == SexStrata, ]%>%
                                               apply(1, mean),
                                             continuous = continuous,
                                           cov.set = c('age'))
                      # 2) aggregation of PRS scores by PCA weights
                      w <- abs(irlba(PRS_X[dfPRS_sample$sex == SexStrata, ],
                                     nv = 1)$v)
                      r2.PCA <- partial.R2(df = dfPRS_sample[dfPRS_sample$sex == SexStrata,],
                                             ypred = (PRS_X %*% w)[dfPRS_sample$sex == SexStrata],
                                             continuous = continuous,
                                           cov.set = c('age'))
                      # 3) the best single PRS score
                      r2.best <- partial.R2(df = dfPRS_sample[dfPRS_sample$sex == SexStrata,],
                                              ypred = PRS_X[dfPRS_sample$sex == SexStrata,
                                                            bestSex$Value[which(bestSex$sex == SexStrata)]],
                                              continuous = continuous,
                                            cov.set = c('age'))
                      # 4) the unsupervised ensemblePRS
                      s.UNSemblePRS <- UNSemblePRS(PRS_X[dfPRS_sample$sex == SexStrata, ])
                      r2.UNSemblePRS <- partial.R2(df = dfPRS_sample[dfPRS_sample$sex == SexStrata, ],
                                                     ypred = s.UNSemblePRS,
                                                     continuous = continuous,
                                                   cov.set = c('age'))
                      data.frame(context = SexStrata,
                                 best = r2.best,
                                 avg = r2.avg,
                                 PCA = r2.PCA,
                                 UNsemblePRS = r2.UNSemblePRS)
                    }) %>% 
                      do.call(rbind, .)
                    
                    rbind(r2.estAncestry, r2.estAge, r2.estSex)
                  }
# stop the parallel
stopCluster(cluster)
r2.estPheno <- data.frame(pheno = pheno, r2.est)
file.name <- paste0('~/', pheno, 
                    '_N', Nsub, 
                    '_niters', niters,
                    '.RData')
save(r2.estPheno, file = file.name)
