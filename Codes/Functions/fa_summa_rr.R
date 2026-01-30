##' Function fa.summa.rr
##'
##' @title Extracting outputs from factor analytic mixed models
##' 
##' @description
##' This function compiles the outputs of factor analytic mixed models fitted using 
##' [asreml::asreml()], rotates the factor loadings and scores using SVD, and builds 
##' the genetic variance-covariance and pairwise genetic correlations matrices.
##' 
##' @param model An object of class `asreml` containing the results of the fitted 
##' factor analytic mixed model.
##' @param env A string that represents the name of the "environment" factor of the model
##' @param gen A string that represents the name of the "genotype" factor of the model
##' 
##' @return The function returns a list with: 
##' \itemize{
##' \item \code{rot.loads}: a matrix with the rotated loadings 
##' \item \code{Gvcov}: the full genetic variance-covariance matrix 
##' \item \code{Gcor}: the matrix of genetic correlations between environments 
##' \item \code{rot.scores}: a matrix with the rotated scores.
##' }
##' 
##' @author Saulo Chaves
##'

fa.summa.rr = function(model, env, geno){
  
  require(asreml)
  
  data = as.data.frame(model$mf)
  temp = attr(attr(model$formulae$random, 'factors'), 'dimnames')[[2]]
  facall = temp[grepl("\\brr\\b", temp) & grepl(geno, temp)]
  envcall = strsplit(facall, split = ':')[[1]][1]
  genocall = strsplit(facall, split = ':')[[1]][2]
  if(!grepl(env, envcall) | !grepl(geno, genocall)){
    envcall = strsplit(facall, split = ':')[[1]][2]
    genocall = strsplit(facall, split = ':')[[1]][1]
  }
  
  num.env = length(model$G.param[[facall]][[envcall]]$levels[-grep("Comp",model$G.param[[facall]][[envcall]]$levels)])
  num.gen = length(model$G.param[[facall]][[genocall]]$levels)
  name.env = model$G.param[[facall]][[envcall]]$levels[-grep("Comp",model$G.param[[facall]][[envcall]]$levels)]
  name.gen = model$G.param[[facall]][[genocall]]$levels
  summa = summary(model)
  
  load = summa$varcomp[grep(paste0("(?=.*", geno,")(?=.*fa\\d+)"), rownames(summa$varcomp), perl = TRUE),]
  var = summa$varcomp[grep(paste0("(?=.*", geno,")(?=.*\\benv:vm\\b)(?!.*\\bvar\\b)(?!.*\\bide\\b)"), rownames(summa$varcomp), perl = TRUE), 1]
  
  
 
  
  # Loadings' rotation via SVD
  load$fa = regmatches(rownames(load), regexpr("fa\\d+", rownames(load)))
  mat.loadings = do.call(cbind, lapply(split(load[,c('fa','component')], f = load$fa), 
                                       function(x) x[,'component']))
  rownames(mat.loadings) = name.env
  
  svd.lambda = svd(mat.loadings)
  D = diag(svd.lambda$d^2, nrow = length(svd.lambda$d))
  if(sum(svd.lambda$u[,1] < 0)/nrow(svd.lambda$u) > 0.5){
    svd.lambda$u = -1 * svd.lambda$u
    svd.lambda$v = -1 * svd.lambda$v
  }
  
  mat.loadings.star = svd.lambda$u
  colnames(mat.loadings.star) = colnames(mat.loadings)
  rownames(mat.loadings.star) = rownames(mat.loadings)
  
  # Scores rotation via SVD
  scores = data.frame(summary(model, coef = T)$coef.random)[
    grep(paste0("(?=.*", geno,")(?=.*Comp)"), rownames(data.frame(summary(model, coef = T)$coef.random)), perl = TRUE),
  ]
  scores$fa = regmatches(rownames(scores), regexpr("Comp\\d+", rownames(scores)))
  scor.vec = c(do.call(rbind, lapply(split(scores, f = scores$fa), function(x) x[,'solution'])))
  scor.vec.star = kronecker(sqrt(D) %*% t(svd.lambda$v), diag(num.gen)) %*% scor.vec
  scor.mat.star = matrix(scor.vec.star, nrow = num.gen, ncol = length(unique(scores$fa)),
                         dimnames = list(name.gen, unique(load$fa)))
  
  
  # Full variance-covariance genetic matrix and genetic correlation matrix
  Gvcov = mat.loadings.star %*% tcrossprod(D, mat.loadings.star) + diag(var)
  Gcor = cov2cor(Gvcov)
  # 
  # # Percentage of explained variance
  # expvar = sum(diag(mat.loadings.star %*% tcrossprod(D, mat.loadings.star)))/
  #   sum(diag(Gvcov))
  # expvar.j = matrix(nrow = ncol(mat.loadings), ncol = num.env, 
  #                   dimnames = list(colnames(mat.loadings), name.env))
  # for (i in 1:nrow(expvar.j)) {
  #   for (j in 1:ncol(expvar.j)) {
  #     expvar.j[i,j] = 100 * mat.loadings.star[j,i]^2 * diag(D)[i]/
  #       (sum(mat.loadings.star[j,]^2 * diag(D)) + var[j])
  #   }
  # }
  # 
  # # Average semivariances ratio
  # lambdacross = mat.loadings.star %*% tcrossprod(D, mat.loadings.star)
  # fullsv = list()
  # for (i in name.env) {
  #   semivar = matrix(nrow = num.env, ncol = 3,
  #                    dimnames = list(name.env,
  #                                    c('i', 'j', 'semivar')))
  #   for (j in name.env) {
  #     semivar[j,] = c(i, j, .5 * (lambdacross[i,i] + lambdacross[j,j]) -
  #                       lambdacross[i,j])
  #   }
  #   fullsv[[i]] = semivar
  #   rm(semivar)
  # }
  # fullsv = data.frame(do.call(rbind, fullsv))
  # fullsv$semivar = as.numeric(fullsv$semivar)
  # fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[,-1]
  # colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  # lambda_ASV = 2/(nrow(fullsv) * (nrow(fullsv)- 1)) * sum(fullsv[upper.tri(fullsv)])
  # rm(fullsv)
  # fullsv = list()
  # for (i in name.env) {
  #   semivar = matrix(nrow = num.env, ncol = 3,
  #                    dimnames = list(name.env,
  #                                    c('i', 'j', 'semivar')))
  #   for (j in name.env) {
  #     semivar[j,] = c(i,j,
  #                     0.5*(Gvcov[i,i] + Gvcov[j,j]) - Gvcov[i,j])
  #   }
  #   fullsv[[i]] = semivar
  # }
  # fullsv = data.frame(do.call(rbind, fullsv))
  # fullsv$semivar = as.numeric(fullsv$semivar)
  # fullsv = reshape(data = fullsv, timevar = 'i', idvar = 'j', direction = 'wide')[,-1]
  # colnames(fullsv) = sub('semivar.', '', colnames(fullsv))
  # G_ASV = 2/(nrow(fullsv) * (nrow(fullsv)- 1)) * sum(fullsv[upper.tri(fullsv)])
  # ASVR = lambda_ASV/G_ASV
  
  results = list('rot.loads' = mat.loadings.star, 
                 'Gvcov' = Gvcov, 
                 'Gcor' = Gcor, 
                 "rot.scores" = scor.mat.star
  )
  # Non-additive variance
  if(any(grepl("ide", temp))){
    nav = summa$varcomp[grep(paste0("(?=.*", geno,")(?=.*\\bide\\b)"), rownames(summa$varcomp), perl = TRUE),]
    results$nav = nav
  }
  return(results)
}

