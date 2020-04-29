#' Random Forest function for genomic prediction.
#'
#' RF_func function will calculate prediction value for input training genotype and phenotype dataset.
#' @param train_pheno the given phenotype data corresponding to genotypes (train_geno)
#' @param train_geno the given genotype data for model training
#' @param test_geno the given genotype data for testing
#' @return A N by 1 vector containing prediction values from test_geno
#' @importFrom stats predict
#' @import randomForest
#' @examples
#'   geno <- matrix(sample(x = 0:2, size = 50000*2000, prob = c(0.45, 0.1, 0.45), replace = TRUE), nrow = 2000)
#'   pheno <- rnorm(2000, mean = 45, sd = 6)
#'   rownames(geno) <- names(pheno) <- paste("Sample", 1:2000)
#'   colnames(geno) <- paste("SNP", 1:50000)
#'   train_samples <- sample(1:2000, 500)
#'   train_geno <- geno[train_samples,]
#'   train_pheno <- pheno[train_samples]
#'   val_geno <- geno[-train_samples,]
#'   val_pheno <- pheno[-train_samples]
#'   pred_pheno <- RF_func(train_pheno, train_geno, val_geno)
#' @export
RF_func <- function(train_pheno, train_geno, test_geno){
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_RF = randomForest::randomForest(x = train_geno, y = train_pheno, verbose=FALSE, mtry=dim(train_geno)[1], ntrees=100, metric="RMSE")
  pred_pheno <- predict(train_RF, newdata = test_geno)
  return_value = list("predicted"=pred_pheno, "model"=train_RF)
  return(return_value)
}
