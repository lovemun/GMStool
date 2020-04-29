#' RRBLUP function for genomic prediction
#'
#' RRb_func function will calculate prediction value for input training genotype and phenotype dataset.
#' @param train_pheno the given phenotype data corresponding to genotypes (train_geno)
#' @param train_geno the given genotype data for model training
#' @param test_geno the given genotype data for testing
#' @return A N by 1 vector containing prediction values from test_geno
#' @importFrom rrBLUP mixed.solve
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
#'   pred_pheno <- RRb_func(train_pheno, train_geno, val_geno)
#' @export
RRb_func = function (train_pheno, train_geno, test_geno) {
  train_geno = (train_geno + 1) / 3
  test_geno = (test_geno + 1) / 3
  train_BLUP <- rrBLUP::mixed.solve(y = train_pheno, Z = train_geno, K = NULL, SE = FALSE, return.Hinv = FALSE)
  train_e = as.matrix(train_BLUP$u)
  pheno_valid = test_geno %*% train_e
  pred_pheno <- pheno_valid + c(train_BLUP$beta)
  return_value = list("predicted"=pred_pheno, "model"=train_BLUP)
  return(return_value)
}
