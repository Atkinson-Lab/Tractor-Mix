#' Run Tractor way of mixed effect model
#'
#' @param fixed formula of your design, this typically include PCs, ages, etc.
#' @param data a data frame contains variable and phenotype in the model.
#' @param kin a positive semi-definite matrix, calculated from PC-Relate from GENESIS.
#' @param id a column in the data frame data, indicating the id of samples.
#' @param family a description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function.
#' @param infiles a vector of input files in order (e.g. anc0.hapcount.txt, anc0.dosage.txt, anc1.dosage.txt); currently only support text file inputs.
#' @param method method of fitting the generalized linear mixed model. Either "REML" or "ML" (default = "REML").
#' @param method.optim optimization method of fitting the generalized linear mixed model. Either "AI", "Brent" or "Nelder-Mead" (default = "AI").
#' @param maxiter a positive integer specifying the maximum number of iterations when fitting the generalized linear mixed model (default = 500).
#' @param tola positive number specifying tolerance, the difference threshold for parameter estimates below which iterations should be stopped. Also the threshold for determining monomorphism. If a SNP has value range less than the tolerance, it will be considered monomorphic and its association test p-value will be NA (default = 1e-5).
#' @param taumin the lower bound of search space for the variance component parameter \tauτ (default = 1e-5), used when method.optim = "Brent".
#' @param taumax the upper bound of search space for the variance component parameter \tauτ (default = 1e5), used when method.optim = "Brent".
#' @param tauregion the number of search intervals for the REML or ML estimate of the variance component parameter \tauτ (default = 10), used when method.optim = "Brent".
#' @param missing.method method of handling missing genotypes. Either "impute2mean" or "omit" (default = "impute2mean").
#' @param verbose a logical switch for printing a progress bar and detailed information (parameter estimates in each iteration) for testing and debugging purpose (default = FALSE).
#' @param ... additional parameters than can be passed to glm
#' @return A list of summary statistics, based on users' input.
#' @export
TractorMix.wald <- function(fixed, data = parent.frame(), kin = NULL, id, family = binomial(link = "logit"), infiles, method = "REML", method.optim = "AI", maxiter = 500, tol = 1e-5, taumin = 1e-5, taumax = 1e5, tauregion = 10, missing.method = "impute2mean", verbose = FALSE, ...) {
  # make sure infiles can take multiple input
  # the first line is header, will be used for sanity check
  # start for the 2nd line

  ################# step 1: read multiple files, and set up an iterator ##############
  # 1. use readLines to read all files
  # 2. save them in a list
  # 3. sanity check if the length is the same

  files = lapply(infiles,function(file){readLines(file)})
  filesName = lapply(infiles, function(file){
    file = unlist(lapply(strsplit(file, "/"), function(vec){vec[length(vec)]}))
    fileSplit = unlist(strsplit(file, ".", fixed = T))
    return(paste(fileSplit[(length(fileSplit)-2):(length(fileSplit)-1)], collapse = "_"))
  })
  names(files) = filesName
  nSNP = length(files[[1]]) - 1


  ################# step 2: setup a result table ##############
  # 1. a list of empty dataframes, each one correspond to each ancestry in file name
  # 2. columns include: CHR, POS, ID, REF, ALT, BETA, SE, P
  resDF = replicate(n = length(files),expr = {setNames(data.frame(matrix(data = NA, nrow = nSNP, ncol = 9)),
                                                       c("CHR", "POS", "ID", "REF", "ALT", "BETA", "SE", "P", "CONVERGED"))},
                    simplify = F)
  names(resDF) = filesName


  ################# step 3: take care of the phenotype and covariates ##############
  # 1. expect an input as the GMMAT example
  # 2. allow users to pass along R formula
  fit0 <- do.call("glm", list(formula = as.formula(deparse(fixed)), data = data, family = family, ...))


  ################# step 4: parse the infiles ##############
  # 1. extract the count information
  # 2. concatenate with the main table
  # 3. fit into glm
  # 4. run wald.test
  for (i in 1:nSNP){
    G = data.frame(lapply(files, function(file){
      fileSplit = unlist(strsplit(file[i+1], "\t"))
      as.integer(fileSplit[6:length(fileSplit)])
    }))
    dataG = cbind(data, G)

    fit0 <- do.call("glm", list(formula = as.formula(paste(deparse(fixed),
                                                           " + ",
                                                           paste(unlist(filesName), collapse = " + "))),
                                data = dataG, family = family, ...))

    group.id = rep(1, nrow(dataG))
    tmpkins = list(kins1 = kin)
    fit <- try(glmmkin.fit(fit0, tmpkins, NULL, group.id, method = method, method.optim = method.optim, maxiter = maxiter, tol = tol, taumin = taumin, taumax = taumax, tauregion = tauregion, verbose = verbose))
    for (j in unlist(filesName)){
      resDF[[j]][i, c("CHR", "POS", "ID", "REF", "ALT")] =  unlist(strsplit(files[[1]][i+1], "\t"))[1:5]
    }


    if(class(fit) != "try-error") {
      # calculate beta, se, and pval for each term
      BETAs = unlist(lapply(seq_along(resDF), function(resdfidx){fit$coefficients[names(resDF)[[resdfidx]]]}))
      SEs = unlist(lapply(seq_along(resDF), function(resdfidx){sqrt(diag(fit$cov)[names(resDF)[[resdfidx]]])}))
      PVALs <- pchisq((BETAs/SEs)^2, 1, lower.tail=F)

      # assign them to the result table
      for (j in unlist(filesName)){
        resDF[[j]][i, "BETA"] = BETAs[j]
        resDF[[j]][i, "SE"] = SEs[j]
        resDF[[j]][i, "P"] = PVALs[j]
        resDF[[j]][i, "CONVERGED"] = fit$converged
      }

    }
  }

  return(resDF)
}


glmmkin.fit <- utils::getFromNamespace("glmmkin.fit", "GMMAT")

