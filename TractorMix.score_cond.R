#' Tractor score test
#'
#' @param obj The output form glmm.kin from GMMAT package
#' @param infiles_geno a vector of input file names in order (e.g. anc0.dosage.txt, anc1.dosage.txt)
#' @param infiles_la a vector of input file names in order (e.g. anc0.hapcount.txt)
#' @param outfiles file name for the summary statistics
#' @param AC_threshold If a variant has *all* ancestry-specific allele counts greater than this value, Tractor-Mix will run a full model; otherwise, Tractor-Mix will drop ancestries with low AC, and only use ancestries with high AC for calculation. 
#' @return A list of summary statistics, based on users' input.
#' @export
#' 

# version: 0.0.1
# date of edit: 02/14/2025
script_version <- "0.0.1"
message(paste("TractorMix.score_cond Script Version :", script_version), "\n")

TractorMix.score_cond <- function(obj, infiles_geno, infiles_la, AC_threshold = 50, outfiles,  n_core = 4, chunk_size = 2048) {
  
  ## install required packages
  required_packages <- c("Matrix", "data.table", "doParallel", "foreach", "dplyr", "abind")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      suppressMessages(library(pkg, character.only = TRUE))
    }
  }
  
  # determine if n_core is correct
  if(n_core>detectCores()) {
    stop(paste0("Error: n_core must be less than ", detectCores(),"!"))
  }
  if (n_core > 1) {
    n_core = n_core - 1
  } else if (n_core <= 0) {
    stop(paste0("Error: --nthreads must be 1 or more!"))
  }
  
  # number of ancestries
  n_anc = length(infiles_geno)
  
  # filter for certain variants 
  iffilter = NA 
  
  # type of test
  if (is.na(AC_threshold)){stop("AC_threshold must be specified")}
  if (AC_threshold %% 1 != 0){stop("AC_threshold must be an integer")}
  
  if (any(grepl("gaussian", as.character(obj$call)))){
    message("Run Tractor-Mix on continuous phenotype")
  } else if (any(grepl("binomial", as.character(obj$call)))){
    message("Run Tractor-Mix on binary phenotype")
  } else {
    stop("You must specify `family = gaussian()` or `family = binomial()` in `glmmkin()`")
  }
  
  
  if (length(infiles_geno) != length(infiles_la) + 1) {
    stop("Conditional analysis requires n genotypes and n - 1 local ancestry terms")
  }
  

  # use wc -l to find the number of SNPs in each file
  # ASSUME: dosage files are equal length
  if (!endsWith(infiles_geno[1],".gz")){
    wc_op <- paste("wc -l", infiles_geno[1]) %>%
      system(intern = T) %>%
      trimws("left") %>% 
      strsplit(" ")
  } else {
    wc_op <- paste("zcat <", infiles_geno[1], "| wc -l") %>%
      system(intern = T) %>%
      trimws("left") %>% 
      strsplit(" ")
  }
  
  
  
  
  n_SNPs = wc_op[[1]][1] %>% as.numeric() -1
  
  
  # setup a result table 
  # 2-way: CHR, POS, ID, REF, ALT, Chi2, P, Eff_anc0, Eff_anc1, Pval_anc0, Pval_anc1, AC_count_anc0, AC_count_anc1, include_anc0, include_anc1
  
  resDF = setNames(data.frame(matrix(data = NA, nrow = 0, ncol = (7 + 5 * n_anc))),
                   c("CHR", "POS", "ID", "REF", "ALT", "Chi2", "P", 
                     paste0("Eff_anc", 0:(n_anc-1)), 
                     paste0("SE_anc", 0:(n_anc-1)), 
                     paste0("Pval_anc", 0:(n_anc-1)),
                     paste0("AC_count", 0:(n_anc-1)),
                     paste0("include_anc", 0:(n_anc-1))))
  
  write.table(resDF, outfiles,  quote = F, row.names = F, sep = "\t")
  
  
  # retrieve file names, and check if consistent
  files_colnames = lapply(infiles_geno, function(infile){colnames(fread(input = infile, sep = "\t", skip = 0, nrows = 1))})

  
  if (all(sapply(files_colnames, function(x) identical(x, files_colnames[[1]])))){
    sample_id = files_colnames[[1]][6:length(files_colnames[[1]])]
  }else{
    stop("Sample order in dosage files do not match!")
  }
  # ASSUME: obj$id_include is a subset of samples in dosage files
  sample2Use = match(obj$id_include, sample_id)
  if (any(is.na(sample2Use))){
    stop("The null model contains samples that are not in the dosage files; 
         remove unmatched samples before fitting the null model")
  }
  
  # define column class
  colclasses = c(rep("character", 5), rep("integer", length(sample_id)))
  
  if(!is.null(obj$P)){
    ################################################################################################
    ################################### Full GRM ###################################################
    ################################################################################################
    X = obj$X
    P = obj$P
    cl <- makeCluster(n_core)
    registerDoParallel(cl)
    pb <- txtProgressBar(min = 0, max = ceiling(n_SNPs / chunk_size), style = 3)
    
    chunk_idx = 0
    while (chunk_idx < ceiling(n_SNPs / chunk_size)){
      files_geno = lapply(infiles_geno, 
                          function(infile){fread(input = infile, 
                                                 sep = "\t", 
                                                 skip = 1 + chunk_idx * chunk_size, 
                                                 nrows = chunk_size, 
                                                 colClasses = colclasses)})
      
      files_la = lapply(infiles_la, 
                        function(infile){fread(input = infile, 
                                               sep = "\t", 
                                               skip = 1 + chunk_idx * chunk_size, 
                                               nrows = chunk_size, 
                                               colClasses = colclasses)})
      
      
      meta_info = files_geno[[1]][,1:5]
      geno_info = abind(lapply(files_geno, function(file){ unname(as.matrix(file[, 6:(length(sample_id) + 5)])[,sample2Use]) }), along = 3)
      la_info = abind(lapply(files_la, function(file){ unname(as.matrix(file[, 6:(length(sample_id) + 5)])[,sample2Use]) }), along = 3)
      
      n_SNPs_chunk = nrow(files_geno[[1]])
      sumstats_chunk = foreach(i = 1:n_SNPs_chunk, .combine = rbind, .inorder = TRUE) %dopar% {
        # parse genotype
        anc_eff = rep(NA, n_anc)
        anc_se = rep(NA, n_anc)
        anc_pval = rep(NA, n_anc)
        
        # parse genotype, filter based on AC 
        G_ = geno_info[i,,]
        LA = la_info[i,,]
        AC = colSums(G_)
        filter_mask = AC > AC_threshold
        
        if (!any(filter_mask)){
          res = c(NA, NA, rep(NA, n_anc * 3), AC, filter_mask)
          names(res) = NULL
          return(res)
        }
        
        iffilter = !all(filter_mask)
        
        # G is the genotype to use, regardless of whether this variants should drop some ancestries
        G = as.matrix(G_[,filter_mask])
        
        Score_G = t(G) %*% obj$scaled.residuals
        Score_LA = t(LA) %*% obj$scaled.residuals
        
        E_ScoreG_cond_LA = t(G) %*% P %*% LA %*% solve(t(LA) %*% P %*%  LA) %*% Score_LA
        Var_ScoreG_cond_LA = t(G) %*% P %*% G - t(G) %*% P %*% LA %*% solve(t(LA) %*% P %*% LA) %*% t(LA) %*% P %*% G
        Score_centered = Score_G - E_ScoreG_cond_LA
        
        if (all(diag(Var_ScoreG_cond_LA) > 0)){
          if (!iffilter){
            joint_chi2 = try(t(Score_centered) %*%  solve(Var_ScoreG_cond_LA) %*% Score_centered, silent=TRUE)
            joint_pval = pchisq(as.numeric(joint_chi2), df = n_anc, lower.tail = F)
            anc_eff = t(Score_centered) %*% solve(Var_ScoreG_cond_LA)
            anc_se = sqrt(diag(solve(Var_ScoreG_cond_LA)))
            anc_pval = sapply((anc_eff/anc_se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
          } else {
            joint_chi2 = NA
            joint_pval = NA
            anc_eff[filter_mask] = t(Score_centered) %*% solve(Var_ScoreG_cond_LA)
            anc_se[filter_mask] = sqrt(diag(solve(Var_ScoreG_cond_LA)))
            anc_pval = sapply((anc_eff/anc_se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
          }
          res = c(round(as.numeric(joint_chi2), 5),
                  signif(as.numeric(joint_pval), 5),
                  round(as.numeric(anc_eff), 5),
                  round(as.numeric(anc_se), 5),
                  signif(as.numeric(anc_pval), 5),
                  AC,
                  filter_mask)
        }else {
          res = c(NA, NA, rep(NA, n_anc * 3), AC, filter_mask)
        }
        names(res) = NULL
        return(res)
        
      }
      write.table(cbind(meta_info, sumstats_chunk), outfiles,  quote = F, row.names = F, sep = "\t", append = T, col.names = F)
      chunk_idx = chunk_idx + 1
      setTxtProgressBar(pb, chunk_idx) 
    }
    close(pb)
    stopCluster(cl)
    
    
  } else {
    ################################################################################################
    ################################### Sparse GRM ###################################################
    ################################################################################################
    
    obj$Sigma_iX = Matrix(obj$Sigma_iX, sparse = TRUE)  # this is solve(Sigma) %*% X
    obj$Sigma_i = Matrix(obj$Sigma_i, sparse = TRUE)  # this is solve(Sigma)
    obj$cov = Matrix(obj$cov, sparse = TRUE)
    
    
    cl <- makeCluster(n_core)
    registerDoParallel(cl)
    pb <- txtProgressBar(min = 0, max = ceiling(n_SNPs / chunk_size), style = 3)
    
    chunk_idx = 0
    while (chunk_idx < ceiling(n_SNPs / chunk_size)){
      files_geno = lapply(infiles_geno, 
                          function(infile){fread(input = infile, 
                                                 sep = "\t", 
                                                 skip = 1 + chunk_idx * chunk_size, 
                                                 nrows = chunk_size, 
                                                 colClasses = colclasses)})
      
      files_la = lapply(infiles_la, 
                        function(infile){fread(input = infile, 
                                               sep = "\t", 
                                               skip = 1 + chunk_idx * chunk_size, 
                                               nrows = chunk_size, 
                                               colClasses = colclasses)})
      
      
      meta_info = files_geno[[1]][,1:5]
      geno_info = abind(lapply(files_geno, function(file){ unname(as.matrix(file[, 6:(length(sample_id) + 5)])[,sample2Use]) }), along = 3)
      la_info = abind(lapply(files_la, function(file){ unname(as.matrix(file[, 6:(length(sample_id) + 5)])[,sample2Use]) }), along = 3)
      
      n_SNPs_chunk = nrow(files_geno[[1]])
      sumstats_chunk = foreach(i = 1:n_SNPs_chunk, .combine = rbind, .inorder = TRUE, .packages=c("Matrix")) %dopar% {
        # parse genotype
        anc_eff = rep(NA, n_anc)
        anc_se = rep(NA, n_anc)
        anc_pval = rep(NA, n_anc)
        
        # parse genotype, filter based on AC 
        G_ = geno_info[i,,]
        LA = la_info[i,,]
        AC = colSums(G_)
        filter_mask = AC > AC_threshold
        
        if (!any(filter_mask)){
          res = c(NA, NA, rep(NA, n_anc * 3), AC, filter_mask)
          names(res) = NULL
          return(res)
        }
        
        iffilter = !all(filter_mask)
        
        # G is the genotype to use, regardless of whether this variants should drop some ancestries
        G = as.matrix(G_[,filter_mask])
        
        Score_G = t(G) %*% obj$scaled.residuals
        Score_LA = t(LA) %*% obj$scaled.residuals
        
        
        Sigma_iG <- obj$Sigma_i %*% G        
        Sigma_iLA <- obj$Sigma_i %*% LA    
        XSigma_iG = t(obj$Sigma_iX) %*% G
        XSigma_iLA = t(obj$Sigma_iX) %*% LA
        
        G_P_LA <- t(G) %*% Sigma_iLA - t(XSigma_iG) %*% obj$cov %*% XSigma_iLA
        LA_P_LA <- t(LA) %*% Sigma_iLA - t(XSigma_iLA) %*% obj$cov %*% XSigma_iLA
        G_P_G <- t(G) %*% Sigma_iG - t(XSigma_iG) %*% obj$cov %*% XSigma_iG
        E_ScoreG_cond_LA <- G_P_LA %*% solve(LA_P_LA) %*% Score_LA
        Var_ScoreG_cond_LA <- G_P_G - G_P_LA %*% solve(LA_P_LA) %*% t(G_P_LA)
        
        # E_ScoreG_cond_LA = t(G) %*% P %*% LA %*% solve(t(LA) %*% P %*%  LA) %*% Score_LA
        # Var_ScoreG_cond_LA = t(G) %*% P %*% G - t(G) %*% P %*% LA %*% solve(t(LA) %*% P %*% LA) %*% t(LA) %*% P %*% G
        Score_centered = Score_G - E_ScoreG_cond_LA
        
        if (all(diag(Var_ScoreG_cond_LA) > 0)){
          if (!iffilter){
            joint_chi2 = try(t(Score_centered) %*%  solve(Var_ScoreG_cond_LA) %*% Score_centered, silent=TRUE)
            joint_pval = pchisq(as.numeric(joint_chi2), df = n_anc, lower.tail = F)
            anc_eff = t(Score_centered) %*% solve(Var_ScoreG_cond_LA)
            anc_se = sqrt(diag(solve(Var_ScoreG_cond_LA)))
            anc_pval = sapply((anc_eff/anc_se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
          } else {
            joint_chi2 = NA
            joint_pval = NA
            anc_eff[filter_mask] = t(Score_centered) %*% solve(Var_ScoreG_cond_LA)
            anc_se[filter_mask] = sqrt(diag(solve(Var_ScoreG_cond_LA)))
            anc_pval = sapply((anc_eff/anc_se)^2, function(teststats){pchisq(as.numeric(teststats), df = 1, lower.tail = F)})
          }
          res = c(round(as.numeric(joint_chi2), 5),
                  signif(as.numeric(joint_pval), 5),
                  round(as.numeric(anc_eff), 5),
                  round(as.numeric(anc_se), 5),
                  signif(as.numeric(anc_pval), 5),
                  AC,
                  filter_mask)
        }else {
          res = c(NA, NA, rep(NA, n_anc * 3), AC, filter_mask)
          
        }
        names(res) = NULL
        return(res)
        
      }
      write.table(cbind(meta_info, sumstats_chunk), outfiles,  quote = F, row.names = F, sep = "\t", append = T, col.names = F)
      chunk_idx = chunk_idx + 1
      setTxtProgressBar(pb, chunk_idx) 
    }
    close(pb)
    stopCluster(cl)
  }
  
  

    
}
  
  
  
  
