# GMStool R package

GWAS-based Marker Selection Tool for Genomic Prediction from Genomic Data. (R >= 3.6.0 is required)

This is GMStool R-package github page. If you are linux or commone line user, please refer to the following github page.

https://github.com/JaeYoonKim72/GMStool

## Requirement libraries
```{r}
require(c("caret", "randomForest", "rrBLUP"))
```

## Installation
```{r}
> devtools::install_github("lovemun/GMStool")
```

## Example
```{r}
> library(GMStool)
> genofile = system.file("data/example_genotype.rds", package = "GMStool")
> phenofile = system.file("data/example_phenotype.rds", package = "GMStool")
> gwasfile = system.file("data/example_gwas.rds", package = "GMStool")
> infofile = system.file("data/Information.txt", package = "GMStool")
> cv = 5; ini_snp = 5; mm = "RRblup"; sel_snps = 1
> load_data = load_files(genofile, phenofile, gwasfile, infofile)

## MAF filtering (optional)
> MAF_QC = geno_QC(load_data, maf_cutoff = 0.05, miss_cutoff = 0.01)

## Cross-validation
> cv_samples = sample(1:cv, nrow(MAF_QC$genotype), replace = TRUE)
> N = ncol(MAF_QC$genotype)

## Marker selection for each cross-validation samples with multi-threads
 # In parallel computing environment, use foreach library for multithreading
> library(tidyverse)
> library(foreach)
> library(iterators)
> results = foreach(j=1:cv %dopar% GMS_main(ini_snps_bk = ini_snp, init_selsnp = sel_snps, j=j, cv_samples = cv_samples,
                    mm = mm, geno2 = MAF_QC$genotype, phenotype1 = MAF_QC$phenotype, preset_fname = NULL,
                    ix = load_data$ix, allm = TRUE, cv = cv, acc1 = 0.9))
 # If not, use for loop
> results = NULL
> for (j in 1:cv){
    results[[j]] = GMS_main(ini_snps_bk = ini_snp, init_selsnp = sel_snps, j=j, cv_samples = cv_samples,
                    mm = mm, geno2 = MAF_QC$genotype, phenotype1 = MAF_QC$phenotype, preset_fname = NULL,
                    ix = load_data$ix, allm = TRUE, cv = cv, acc1 = 0.9)
  }

## Choose final markers and generate model
> all_train_acc <- NULL; selected_train_acc <- NULL
> CV_results <- list(); nSamGenSum <- list(); tsum <- list(); inisum <- list()
> results = get(x = paste0("results_", mm))
> for (k in 1:(cv-1)){
	all_train_acc <- c(all_train_acc, results[[k]]$all_train_acc)
	selected_train_acc <- c(selected_train_acc, results[[k]]$selected_train_acc)
	CV_results[[paste0("CV-", k, "_", mm)]] <- results[[k]]$CV_results
	nSamGenSum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$nSamGenSum
	tsum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$tsum
	inisum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$inisum
}
