#' Loading input file function for Marker selection
#'
#' load_inputs function will load datasets using only one command.
#' @param genofile the given genotype data file path (Required)
#' @param phenofile the given phenotype data file path (Required)
#' @param mmm the given method for genomic prediction
#' @param gwasfile the given GWAS result data file path. Default is NULL.
#' @param infofile the given SNP information file path. Default is NULL.
#' @return A list of all input files
#' @importFrom data.table fread
#' @importFrom rrBLUP mixed.solve
#' @importFrom caret train varImp
#' @importFrom utils read.table write.csv
#' @export
load_files = function(genofile, phenofile, gwasfile = NULL, infofile = NULL, mmm = "RRblup"){
  if (endsWith(x = genofile, suffix = "rds")){
    geno <- readRDS(genofile)
  } else {
    geno <- data.table::fread(genofile, header=T, check.names=FALSE)
  }
  SNP_ids <- as.vector(unlist(geno[,1]))
  sample_ids <- colnames(geno)[-1]
  geno <- t(geno[,-1])
  colnames(geno) <- SNP_ids
  rownames(geno) <- sample_ids

  if (endsWith(x = phenofile, suffix = "rds")){
    pheno <- readRDS(phenofile)
  } else {
    pheno <- read.table(phenofile, header=T, check.names=F)
  }
  rownames(pheno) <- pheno[,1]
  common_samples <- intersect(rownames(geno), rownames(pheno))
  geno1 <- geno[common_samples,]
  pheno1 <- data.frame(pheno[common_samples,], check.names=F)
  rm(geno); rm(pheno)
  pheno_num <- 2
  cat(paste0("Phenotype: Name of phenotype is ", paste(colnames(pheno1)[pheno_num], '\n', sep="")))

  phenotype <- pheno1[,pheno_num]
  names(phenotype) <- rownames(pheno1)
  pheno_NA <- which(is.na(phenotype))
  if (length(pheno_NA) > 0){
    geno1 <- geno1[-pheno_NA,]
    phenotype <- phenotype[-pheno_NA]
    pheno1 <- data.frame(pheno1[-pheno_NA,], check.names=F)
  }
  cat(paste0("- # of missing phenotype is ", length(pheno_NA)), '\n\n')

  #(3)Load GWAS results or calcuate marker effects
  if (is.null(gwasfile)){
    cat("GWAS: Calculate genomic marker effects.", '\n')
    cat("- It takes a lot of time.", '\n')
    if (mmm == "RRblup" | mmm == "RRblup_RF"){
      BLUP = rrBLUP::mixed.solve(y = phenotype, Z = geno1, K = NULL, SE = FALSE, return.Hinv = FALSE)
      marker_effects = BLUP$u
    } else if (mmm == "RF"){
      ttMod <- caret::train(x = geno1, y = phenotype, method = "rf")
      marker_effects = varImp(ttMod)$importance
    }
    gwas_results <- data.frame(SNP = names(marker_effects), P = marker_effects)
    gwas_results <- gwas_results[order(abs(gwas_results$P), decreasing = TRUE),]
    gwas_results <- gwas_results[gwas_results$SNP == (intersect(gwas_results$SNP, colnames(geno1))),]
    #gwas_results <- gwas_results[which(gwas_results$SNP %in% colnames(geno1)),]
    final_snps <- as.vector(gwas_results$SNP)
    geno1 <- geno1[,final_snps]
    cat(paste0("GWAS: Saving GWAS ordered by marker effect (GWAS_marker_effect_ordered.csv)", '\n\n'), sep="")
    write.csv(gwas_results, file = "GWAS_marker_effect_ordered.csv", quote=F)
    ix <- final_snps
    if (!is.null(infofile)){
      infotable <- fread(infofile, header=T, check.names=FALSE, data.table=FALSE)
    }
  } else {
    snp_cn=1
    chr_cn=2
    pos_cn=3
    pv_cn=4
    if (endsWith(x = gwasfile, suffix = "rds")){
      gwas_results <- readRDS(gwasfile)
    } else {
      gwas_results <- fread(gwasfile, header=T, check.names=FALSE, data.table=FALSE)
    }
    write.csv(gwas_results, file = "GAWS_marker_pvalue.csv", quote=F, row.names=F)
    gwas_results <- gwas_results[order(gwas_results[,pv_cn]),]
    gwas_results <- gwas_results[gwas_results[,snp_cn] == (intersect(gwas_results[,snp_cn], colnames(geno1))),]
    final_snps <- as.vector(gwas_results[,snp_cn])
    geno1 <- geno1[,final_snps]
    cat(paste0("GWAS: Saving GWAS result ordered by P-value (GWAS_marker_pvalue_ordered.csv).", '\n\n'), sep="")
    write.csv(gwas_results, file = "GWAS_marker_pvalue_ordered.csv", quote=F)
    ix <- final_snps
  }

  if (dim(geno1)[2] == 0){  cat(paste0("Error3: Marker names is not matched between  Genotype file and GWAS file."))
    quit(save="no")
  }
  if (is.null(gwasfile)){
    return(list(genotype = geno1, phenotype = phenotype, pheno = pheno1, gwas = gwas_results, info = infotable, ix = ix))
  } else {
    return(list(genotype = geno1, phenotype = phenotype, pheno = data.frame(pheno1), gwas = gwas_results, ix = ix))
  }
}
