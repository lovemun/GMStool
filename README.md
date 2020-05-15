# GMStool R package

GWAS-based Marker Selection Tool for Genomic Prediction from Genomic Data. (R >= 3.6.0 is required)

This is GMStool R-package github page. If you are linux or commone line user, please refer to the following github page.

https://github.com/JaeYoonKim72/GMStool

## Requirement libraries
```{r}
requirements <- c("caret", "randomForest", "rrBLUP")
new.packages <- requirements[!(requirements %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
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
> cv = 5; ini_snp = 5; mm = "RRblup"; sel_snps = 1; preset_fname = NULL
> load_data = load_files(genofile, phenofile, gwasfile, infofile)

## MAF filtering (optional)
> MAF_QC = geno_QC(load_data, maf_cutoff = 0.05, miss_cutoff = 0.01)

## Cross-validation
> cv_samples = sample(1:cv, nrow(MAF_QC$genotype), replace = TRUE)
> N = ncol(MAF_QC$genotype)
> t_time <- Sys.time()

## Marker selection for each cross-validation samples with multi-threads
 # In parallel computing environment, use foreach library for multithreading
> library(tidyverse)
> library(foreach)
> library(iterators)
> results = foreach(j=1:cv %dopar% GMS_main(ini_snps_bk = ini_snp, init_selsnp = sel_snps, j=j, cv_samples = cv_samples,
                    mm = mm, geno2 = MAF_QC$genotype, phenotype1 = MAF_QC$phenotype, preset_fname = NULL,
                    ix = load_data$ix, allm = TRUE, cv = cv, acc1 = 0.9, time_cutoff = 1, delta = 0.001))
 # If not, use for loop
> results = NULL
> for (j in 1:cv){
    results[[j]] = GMS_main(ini_snps_bk = ini_snp, init_selsnp = sel_snps, j=j, cv_samples = cv_samples,
                    mm = mm, geno2 = MAF_QC$genotype, phenotype1 = MAF_QC$phenotype, preset_fname = NULL,
                    ix = load_data$ix, allm = TRUE, cv = cv, acc1 = 0.9, time_cutoff = 1, delta = 0.001)
  }

## Choose final markers and generate model
> all_train_acc <- NULL; selected_train_acc <- NULL
> CV_results <- list(); nSamGenSum <- list(); tsum <- list(); inisum <- list()
> for (k in 1:cv){
	all_train_acc <- c(all_train_acc, results[[k]]$all_train_acc)
	selected_train_acc <- c(selected_train_acc, results[[k]]$selected_train_acc)
	CV_results[[paste0("CV-", k, "_", mm)]] <- results[[k]]$CV_results
	nSamGenSum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$nSamGenSum
	tsum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$tsum
	inisum[[paste0("CV-", k, "_", mm)]] <- results[[k]]$inisum
}

## Save results and summary statistics
> CV_results = CV_results[order(names(CV_results), decreasing=FALSE)]
> nSamGenSum = nSamGenSum[order(names(nSamGenSum), decreasing=FALSE)]
> inisum = inisum[order(names(inisum), decreasing=FALSE)]
> tsum = tsum[order(names(tsum), decreasing=FALSE)]
> marker_vec <- unlist(CV_results)
> marker_table <- table(marker_vec)[order(table(marker_vec), decreasing = TRUE)]
> marker_select <- names(marker_table)
> msummary_list <- list()
> name_tmp = paste("CV-", seq(cv), "_", mm, sep="")
> for (iter in names(CV_results)){
  msummary_list[[iter]] = sapply(marker_select, function(x) as.integer(x %in% CV_results[[iter]]))
}
> marker_summary = data.frame(msummary_list, check.names = FALSE)
> marker_summary = marker_summary[,name_tmp]
> marker_summary = data.frame("Num"=seq(1,dim(marker_summary)[1]), "Marker"=rownames(marker_summary), marker_summary, "Total_score"=apply(marker_summary,1,sum), check.names = FALSE)
> write.table(marker_summary, file=paste("CV_", mm, "_Marker_scores.txt", sep=""), quote=FALSE, sep="\t", row.names = FALSE)
> write(marker_select, file = paste("CV_Marker_lists.txt", sep=""))
> m_all_train_acc = round(mean(all_train_acc),6)
> sd_all_train_acc = round(sd(all_train_acc),6)
> add_row = c("Total", mm, dim(MAF_QC$genotype)[1], dim(MAF_QC$genotype)[1], dim(MAF_QC$genotype)[2], 
	paste(round(mean(all_train_acc), 6), " (",round(sd(all_train_acc),6),")", sep="") , 
	paste(round(mean(selected_train_acc), 6), " (",round(sd(selected_train_acc),6),")", sep=""), length(preset_fname),
	paste(length(unique(marker_vec)), " (", sum(table(marker_vec) == cv), ")", sep=""), "-", 
	format(difftime(Sys.time(), t_time), usetz = TRUE)) 
> names(all_train_acc) = name_tmp; names(selected_train_acc) = name_tmp
> tsummary = data.frame("Corr_all_markers"=all_train_acc[names(CV_results)],
	"Corr_selected_markers"=selected_train_acc[names(CV_results)])
> tsummary = data.frame("CV"=unlist(lapply(strsplit(names(CV_results), "_"), function(x) x[1])),
	"Model"=unlist(lapply(strsplit(names(CV_results), "_"), function(x) x[2])), 
	t(data.frame(nSamGenSum, check.names = FALSE))[names(CV_results),], tsummary, "Pre-selected_markers"=length(preset_fname),
	"Total-Selected_markers"=unlist(lapply(CV_results, length)), t(data.frame(inisum, check.names = FALSE))[names(CV_results),],
	t(data.frame(tsum, check.names = FALSE))[names(CV_results),])
> colnames(tsummary) = c("CV", "Model", "Train_Samples", "Test_Samples", "Used_Markers", "Corr_all_markers","Corr_selected_markers",
	"Pre_selected_markers", "Total_Selected_markers", "Initial_runtime", "Total_runtime")
> tsummary = tsummary[order(tsummary$Model),]
> tsummary = rbind(as.matrix(tsummary), add_row)
> write.table(tsummary, file=paste("CV_", mm, "_Selection_summary.txt", sep=""), quote=FALSE, sep="\t", row.names = FALSE)
