##' Function expvar
##'
##' @title Compute the variance explained by a Factor Analytic Mixed Model
##' 
##' @description
##' This function uses the factor loadings to estimate the explicate power of 
##' factor analytic mixed models fitted using [asreml::asreml()]
##' 
##' @param model An object of class `asreml` containing the results of the fitted 
##' factor analytic mixed model.
##' @param env A string that represents the name of the "environment" factor of the model
##' @param gen A string that represents the name of the "genotype" factor of the model
##' 
##' @return The function returns a data frame with the Akaike Information Criterion, 
##' the overall variance explained by the model, the variance explained by each factor, 
##' and the variance explained of each environment.
##' 
##' @author Saulo Chaves
##'

expvar = function(famod, env, geno){
  env.num = nlevels(famod$mf[,env])
  env.names = levels(famod$mf[,env])
  # geno.num = length(unique(geno.names))
  vc = summary(famod)$varcomp
  fa = unique(na.exclude(str_extract(rownames(vc), "fa[0-9]")))
  L = matrix(vc[grep('fa', rownames(vc)),1], nrow = env.num, ncol = length(fa),
             dimnames = list(env.names, fa))
  Psi = diag(vc[which(grepl(geno, rownames(vc)) & !grepl("rr|fa|ide", rownames(vc))), 1])
  svd.L = svd(L)
  if(sum(svd.L$u[,1] < 0)/nrow(svd.L$u) >= 0.5){
    svd.L$u = -1 * svd.L$u
    svd.L$v = -1 * svd.L$v
  }
  L.star = svd.L$u
  D = diag(x = svd.L$d^2, nrow = length(svd.L$d), ncol = length(svd.L$d))
  G = L.star %*% tcrossprod(D, L.star) + Psi
  expvar = sum(diag(D))/sum(diag(G)) * 100
  expvar.k = diag(D)/sum(diag(G)) * 100
  expvar.j = diag(L.star %*% tcrossprod(D, L.star))/diag(G) * 100
  
  return(data.frame(
    level = c("AIC","Average", paste0("FA", seq_along(expvar.k)), env.names),
    level2 = c("..", "..", rep("k", length(expvar.k)), rep("j",env.num)),
    expvar = c(summary(famod)$aic, expvar, expvar.k, expvar.j),
    mod = paste0("FA", length(expvar.k))
  ))
}
