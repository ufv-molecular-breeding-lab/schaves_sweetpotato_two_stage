source("codes/0_preps.R")
load(file = "saves/resp2rev/stg1_models.RDA")

nfolds = 5
nrepl = 30


# 2S ----------------------------------------------------------------------

# BLUE --------------------------------------------------------------------

## Both pools --------------------------------------------------------------
blue = do.call(rbind, lapply(models, function(x) x$BLUE$pred.df))
V_blue = do.call(bdiag, lapply(models, function(x) x$BLUE$bdVt))
V_blue = V_blue[which(grepl("UGP", blue$geno.cod)), which(grepl("UGP", blue$geno.cod))]
blue = droplevels(blue[which(grepl("UGP", blue$geno.cod)), ])

ginv = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blue$geno.cod),
                          colnames(Gmat) %in% levels(blue$geno.cod)],
                 sparseform = TRUE)$Ginv.sparse

V_blue = V_blue[which(blue$geno.cod %in% attr(ginv, 'rowNames')), which(blue$geno.cod %in% attr(ginv, 'rowNames'))]
blue = droplevels(blue[which(blue$geno.cod %in% attr(ginv, 'rowNames')),])
blue$TG = as.factor(seq_along(blue$geno.cod))
V_blue@Dimnames = list(as.character(blue$TG), as.character(blue$TG))
V_blue_sparse = full2sparse(as.matrix(V_blue))
attr(V_blue_sparse, 'INVERSE') = FALSE
blue$env = as.factor(blue$env)
aux = na.omit(unique(blue[,c("geno.cod",'pool')]))

V_diag = matrix(0, nrow = nrow(V_blue), ncol = ncol(V_blue))
diag(V_diag) = 1/diag(solve(V_blue))
rownames(V_diag) = rownames(V_blue)
colnames(V_diag) = colnames(V_blue)
attr(V_diag, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 7
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = blue
  cvdata[[rept]]$set = NA
  for (id in blue$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blue", "cvdata", "nfolds", "nrepl", "V_diag"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV2 - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_diag),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod)+ vm(TG, V_diag),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV2_DW.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blue", "cvdata", "nfolds", "nrepl", "V_blue_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV2 - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV2_FW.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(blue$geno)) * nrepl
  )[order(runif(length(unique(blue$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = blue
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(blue$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blue", "cvdata", "nfolds", "nrepl",'V_diag'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV1 - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_diag),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod)+ vm(TG, V_diag),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV1_DW.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blue", "cvdata", "nfolds", "nrepl", "V_blue_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV1 - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV1_FW.RDS")


## Pool A --------------------------------------------------------------
V_blueA = V_blue[which(blue$pool == "A"), which(blue$pool == "A")]
blueA = droplevels(blue[which(blue$pool == "A"), ])
ginvA = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blueA$geno.cod), colnames(Gmat) %in% levels(blueA$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_blue_sparse = full2sparse(as.matrix(V_blueA))
V_diagA = matrix(0, nrow = nrow(V_blueA), ncol = ncol(V_blueA))
diag(V_diagA) = 1/diag(solve(V_blueA))
rownames(V_diagA) = rownames(V_blueA)
colnames(V_diagA) = colnames(V_blueA)
attr(V_diagA, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 7
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = blueA
  cvdata[[rept]]$set = NA
  for (id in blueA$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blueA", "cvdata", "nfolds", "nrepl","V_diagA"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV2 (A) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagA),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod)  + vm(TG, V_diagA),
        data = w,
        family = asr_gaussian(dispersion = 0.001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV2_DW_poolA.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blueA", "cvdata", "nfolds", "nrepl", "V_blue_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV2 (A) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV2_FW_poolA.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(blueA$geno)) * nrepl
  )[order(runif(length(unique(blueA$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = blueA
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(blueA$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blueA", "cvdata", "nfolds", "nrepl", "V_diagA"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV1 (A) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagA),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod)  + vm(TG, V_diagA),
        data = w,
        family = asr_gaussian(dispersion = 0.001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV1_DW_poolA.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blueA", "cvdata", "nfolds", "nrepl", "V_blue_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV1 (A) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV1_FW_poolA.RDS")


## Pool B --------------------------------------------------------------
V_blueB = V_blue[which(blue$pool == "B"), which(blue$pool == "B")]
blueB = droplevels(blue[which(blue$pool == "B"), ])
ginvB = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blueB$geno.cod), colnames(Gmat) %in% levels(blueB$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_blue_sparse = full2sparse(as.matrix(V_blueB))

V_diagB = matrix(0, nrow = nrow(V_blueB), ncol = ncol(V_blueB))
diag(V_diagB) = 1/diag(solve(V_blueB))
rownames(V_diagB) = rownames(V_blueB)
colnames(V_diagB) = colnames(V_blueB)
attr(V_diagB, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 77
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = blueB
  cvdata[[rept]]$set = NA
  for (id in blueB$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blueB", "cvdata", "nfolds", "nrepl",'V_diagB'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV2 (B) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagB),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_diagB),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV2_DW_poolB.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blueA", "cvdata", "nfolds", "nrepl", "V_blue_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV2 (B) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV2_FW_poolB.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(77)
sets = split(
  rep(
    1:nfolds, length(unique(blueB$geno)) * nrepl
  )[order(runif(length(unique(blueB$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = blueB
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(blueB$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blueB", "cvdata", "nfolds", "nrepl", "V_diagB"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV1 (B) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagB),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_diagB),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV1_DW_poolB.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blueA", "cvdata", "nfolds", "nrepl", "V_blue_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("BLUE - CV1 (B) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/BLUE_CV1_FW_poolB.RDS")



# dBLUP --------------------------------------------------------------------
## Both pools --------------------------------------------------------------
blup = do.call(rbind, lapply(models, function(x) x$BLUP$pred.df))
V_blup = do.call(bdiag, lapply(models, function(x) x$BLUP$bdVt))
V_blup = V_blup[which(paste(blup$geno.cod, blup$pool, sep = '@') %in%
                        paste(aux$geno.cod, aux$pool, sep = '@')), which(paste(blup$geno.cod, blup$pool, sep = '@') %in%
                                                                           paste(aux$geno.cod, aux$pool, sep = '@'))]
blup = droplevels(blup[which(paste(blup$geno.cod, blup$pool, sep = '@') %in%
                               paste(aux$geno.cod, aux$pool, sep = '@')), ])
V_blup = V_blup[which(paste(blup$geno.cod, blup$env, sep = "@") %in% paste(blue$geno.cod, blue$env, sep = "@")),
                which(paste(blup$geno.cod, blup$env, sep = "@") %in% paste(blue$geno.cod, blue$env, sep = "@"))]
blup = droplevels(blup[which(paste(blup$geno.cod, blup$env, sep = "@") %in% paste(blue$geno.cod, blue$env, sep = "@")),])

ginv = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blup$geno.cod),
                          colnames(Gmat) %in% levels(blup$geno.cod)],
                 sparseform = TRUE)$Ginv.sparse

blup$TG = as.factor(seq_along(blup$geno.cod))
V_blup@Dimnames = list(as.character(blup$TG), as.character(blup$TG))
V_blup_sparse = full2sparse(as.matrix(V_blup))
attr(V_blup_sparse, 'INVERSE') = FALSE
blup$env = as.factor(blup$env)

V_diag = matrix(0, nrow = nrow(V_blup), ncol = ncol(V_blup))
diag(V_diag) = 1/diag(solve(V_blup))
rownames(V_diag) = rownames(V_blup)
colnames(V_diag) = colnames(V_blup)
attr(V_diag, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 7
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = blup
  cvdata[[rept]]$set = NA
  for (id in blup$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blup", "cvdata", "nfolds", "nrepl",'V_diag'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV2 - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_diag),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod)+ vm(TG, V_diag),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV2_DW.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blup", "cvdata", "nfolds", "nrepl", "V_blup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV2 - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV2_FW.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(blup$geno)) * nrepl
  )[order(runif(length(unique(blup$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = blup
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(blup$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blup", "cvdata", "nfolds", "nrepl", "V_diag"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV1 - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_diag),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod)+ vm(TG, V_diag),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV1_DW.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "blup", "cvdata", "nfolds", "nrepl", "V_blup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV1 - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
          corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV1_FW.RDS")


## Pool A --------------------------------------------------------------
V_blupA = V_blup[which(blup$pool == "A"), which(blup$pool == "A")]
blupA = droplevels(blup[which(blup$pool == "A"), ])
ginvA = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blupA$geno.cod), colnames(Gmat) %in% levels(blupA$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_blup_sparse = full2sparse(as.matrix(V_blupA))
V_diagA = matrix(0, nrow = nrow(V_blupA), ncol = ncol(V_blupA))
diag(V_diagA) = 1/diag(solve(V_blupA))
rownames(V_diagA) = rownames(V_blupA)
colnames(V_diagA) = colnames(V_blupA)
attr(V_diagA, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 7
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = blupA
  cvdata[[rept]]$set = NA
  for (id in blupA$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blupA", "cvdata", "nfolds", "nrepl",'V_diagA'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV2 (A) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagA),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod) + vm(TG, V_diagA),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV2_DW_poolA.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blupA", "cvdata", "nfolds", "nrepl", "V_blup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV2 (A) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV2_FW_poolA.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(blupA$geno)) * nrepl
  )[order(runif(length(unique(blupA$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = blupA
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(blupA$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blupA", "cvdata", "nfolds", "nrepl",'V_diagA'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV1 (A) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagA),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod) + vm(TG, V_diagA),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV1_DW_poolA.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "blupA", "cvdata", "nfolds", "nrepl", "V_blup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV1 (A) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
          corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV1_FW_poolA.RDS")


## Pool B --------------------------------------------------------------
V_blupB = V_blup[which(blup$pool == "B"), which(blup$pool == "B")]
blupB = droplevels(blup[which(blup$pool == "B"), ])
ginvB = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blupB$geno.cod), colnames(Gmat) %in% levels(blupB$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_blup_sparse = full2sparse(as.matrix(V_blupB))

V_diagB = matrix(0, nrow = nrow(V_blupB), ncol = ncol(V_blupB))
diag(V_diagB) = 1/diag(solve(V_blupB))
rownames(V_diagB) = rownames(V_blupB)
colnames(V_diagB) = colnames(V_blupB)
attr(V_diagB, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 77
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = blupB
  cvdata[[rept]]$set = NA
  for (id in blupB$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blupB", "cvdata", "nfolds", "nrepl",'V_diagB'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV2 (B) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagB),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_diagB),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV2_DW_poolB.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blupB", "cvdata", "nfolds", "nrepl", "V_blup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV2 (B) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV2_FW_poolB.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(77)
sets = split(
  rep(
    1:nfolds, length(unique(blupB$geno)) * nrepl
  )[order(runif(length(unique(blupB$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = blupB
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(blupB$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blupB", "cvdata", "nfolds", "nrepl","V_diagB"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV1 (B) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagB),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_diagB),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV1_DW_poolB.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "blupB", "cvdata", "nfolds", "nrepl", "V_blup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dBLUP - CV1 (B) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
          corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dBLUP_CV1_FW_poolB.RDS")


# dABLUP --------------------------------------------------------------------
## Both pools --------------------------------------------------------------
ablup = do.call(rbind, lapply(models, function(x) x$ABLUP$pred.df))
V_ablup = do.call(bdiag, lapply(models, function(x) x$ABLUP$bdVt))
V_ablup = V_ablup[which(paste(ablup$geno.cod, ablup$pool, sep = '@') %in%
                          paste(aux$geno.cod, aux$pool, sep = '@')), which(paste(ablup$geno.cod, ablup$pool, sep = '@') %in%
                                                                             paste(aux$geno.cod, aux$pool, sep = '@'))]
ablup = droplevels(ablup[which(paste(ablup$geno.cod, ablup$pool, sep = '@') %in%
                                 paste(aux$geno.cod, aux$pool, sep = '@')), ])
V_ablup = V_ablup[which(paste(ablup$geno.cod, ablup$env, sep = "@") %in% paste(blue$geno.cod, blue$env, sep = "@")),
                  which(paste(ablup$geno.cod, ablup$env, sep = "@") %in% paste(blue$geno.cod, blue$env, sep = "@"))]
ablup = droplevels(ablup[which(paste(ablup$geno.cod, ablup$env, sep = "@") %in% paste(blue$geno.cod, blue$env, sep = "@")),])

ginv = G.inverse(G = Gmat[rownames(Gmat) %in% levels(ablup$geno.cod),
                          colnames(Gmat) %in% levels(ablup$geno.cod)],
                 sparseform = TRUE)$Ginv.sparse

ablup$TG = as.factor(seq_along(ablup$geno.cod))
V_ablup@Dimnames = list(as.character(ablup$TG), as.character(ablup$TG))
V_ablup_sparse = full2sparse(as.matrix(V_ablup))
attr(V_ablup_sparse, 'INVERSE') = FALSE
ablup$env = as.factor(ablup$env)
V_diag = matrix(0, nrow = nrow(V_ablup), ncol = ncol(V_ablup))
diag(V_diag) = 1/diag(solve(V_ablup))
rownames(V_diag) = rownames(V_ablup)
colnames(V_diag) = colnames(V_ablup)
attr(V_diag, 'INVERSE') = FALSE
### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 7
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = ablup
  cvdata[[rept]]$set = NA
  for (id in ablup$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "ablup", "cvdata", "nfolds", "nrepl","V_diag"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV2 - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_diag),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) +
          corh(env):ide(geno.cod) + vm(TG, V_diag),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV2_DW.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "ablup", "cvdata", "nfolds", "nrepl", "V_ablup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV2 - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv)  + 
        corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv)  + 
          corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV2_FW.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(ablup$geno)) * nrepl
  )[order(runif(length(unique(ablup$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = ablup
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(ablup$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "ablup", "cvdata", "nfolds", "nrepl","V_diag"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV1 - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
        corh(env):ide(geno.cod) + vm(TG, V_diag),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) +
          corh(env):ide(geno.cod) + vm(TG, V_diag),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV1_DW.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "ablup", "cvdata", "nfolds", "nrepl", "V_ablup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV1 - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) +
        corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 2):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) +
          corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 50
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginv\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV1_FW.RDS")


## Pool A --------------------------------------------------------------
V_ablupA = V_ablup[which(ablup$pool == "A"), which(ablup$pool == "A")]
ablupA = droplevels(ablup[which(ablup$pool == "A"), ])
ginvA = G.inverse(G = Gmat[rownames(Gmat) %in% levels(ablupA$geno.cod), colnames(Gmat) %in% levels(ablupA$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_ablup_sparse = full2sparse(as.matrix(V_ablupA))

V_diagA = matrix(0, nrow = nrow(V_ablupA), ncol = ncol(V_ablupA))
diag(V_diagA) = 1/diag(solve(V_ablupA))
rownames(V_diagA) = rownames(V_ablupA)
colnames(V_diagA) = colnames(V_ablupA)
attr(V_diagA, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 7
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = ablupA
  cvdata[[rept]]$set = NA
  for (id in ablupA$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "ablupA", "cvdata", "nfolds", "nrepl","V_diagA"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV2 (A) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagA),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA)+
          corh(env):ide(geno.cod) + vm(TG, V_diagA),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV2_DW_poolA.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "ablupA", "cvdata", "nfolds", "nrepl", "V_ablup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV2 (A) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) +
        corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) +
          corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV2_FW_poolA.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(ablupA$geno)) * nrepl
  )[order(runif(length(unique(ablupA$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = ablupA
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(ablupA$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "ablupA", "cvdata", "nfolds", "nrepl",'V_diagA'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV1 (A) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagA),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA)+
          corh(env):ide(geno.cod) + vm(TG, V_diagA),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV1_DW_poolA.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvA", "ablupA", "cvdata", "nfolds", "nrepl", "V_ablup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV1 (A) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) +
        corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + vm(TG, V_ablup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvA\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV1_FW_poolA.RDS")


## Pool B --------------------------------------------------------------
V_ablupB = V_ablup[which(ablup$pool == "B"), which(ablup$pool == "B")]
ablupB = droplevels(ablup[which(ablup$pool == "B"), ])
ginvB = G.inverse(G = Gmat[rownames(Gmat) %in% levels(ablupB$geno.cod), colnames(Gmat) %in% levels(ablupB$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_ablup_sparse = full2sparse(as.matrix(V_ablupB))

V_diagB = matrix(0, nrow = nrow(V_ablupB), ncol = ncol(V_ablupB))
diag(V_diagB) = 1/diag(solve(V_ablupB))
rownames(V_diagB) = rownames(V_ablupB)
colnames(V_diagB) = colnames(V_ablupB)
attr(V_diagB, 'INVERSE') = FALSE

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 77
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = ablupB
  cvdata[[rept]]$set = NA
  for (id in ablupB$geno) {
    cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(1:nfolds,
                                                              size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
                                                              replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

#### Diagonal ---------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "ablupB", "cvdata", "nfolds", "nrepl","V_diagB"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV2 (B) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagB),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB)+
          corh(env):ide(geno.cod) + vm(TG, V_diagB),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV2_DW_poolB.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "ablupB", "cvdata", "nfolds", "nrepl", "V_ablup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV2 (B) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) +
        corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + vm(TG, V_ablup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV2_FW_poolB.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(77)
sets = split(
  rep(
    1:nfolds, length(unique(ablupB$geno)) * nrepl
  )[order(runif(length(unique(ablupB$geno)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = ablupB
  cvdata = merge(cvdata, data.frame(
    geno.cod = unique(ablupB$geno.cod),
    set = x
  ), by = 'geno.cod')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

#### Diagonal ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "ablupB", "cvdata", "nfolds", "nrepl",'V_diagB'))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV1 (B) - Diagonal // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = rytha ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
        corh(env):ide(geno.cod) + vm(TG, V_diagB),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    
    mod = tryCatch({
      asreml(
        fixed = yNA ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB)+
          corh(env):ide(geno.cod) + vm(TG, V_diagB),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        workspace = "2gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV1_DW_poolB.RDS")

#### Full ----------------------------------------------------------------
message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginvB", "ablupB", "cvdata", "nfolds", "nrepl", "V_ablup_sparse"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("dABLUP - CV1 (B) - Full // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = asreml(
      fixed = yNA ~ -1 + env,
      random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) +
        corh(env):ide(geno.cod)+ vm(TG, V_ablup_sparse),
      data = w,
      family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
    )
    params = mod$vparameters.table
    params$Value[grep("TG", params$Component)] = 1
    params$Constraint[grep("TG", params$Component)] = "F"
    
    mod = tryCatch({
      asreml(
        fixed = rytha ~ -1 + env,
        random = ~ rr(env, 1):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB)+
          corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
        data = w,
        family = asr_gaussian(dispersion = 0.0001),
        na.action = na.method(x = 'exclude', y = 'exclude'),
        G.param = params,
        workspace = "3.5gb",
        maxit = 100
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      
      blu = summary(mod, coef = TRUE)$coef.random
      rr = data.frame(marg = blu[which(
        grepl("vm", rownames(blu)) &
          grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      rr = rr |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("rr\\(env, [0-9]\\)_", '', env))
      delta = data.frame(delta = blu[which(
        grepl("vm", rownames(blu)) & grepl("geno.cod", rownames(blu)) &
          !grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1])
      delta = delta |> rownames_to_column('env') |>
        separate('env', into = c("env", 'geno.cod'), sep = ':vm\\(geno.cod, ginvB\\)_') |>
        mutate(env = gsub("env_", '', env))
      predvalu[[j]] = left_join(rr, delta, by = c("env", 'geno.cod')) |> 
        mutate(yhat = marg + delta) |> 
        right_join(droplevels(w[which(is.na(w$yNA)), ]), by = c('geno.cod', 'env')) |> 
        dplyr::select(env, geno.cod, yhat, rytha, set, rept)
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/dABLUP_CV1_FW_poolB.RDS")


# SS ----------------------------------------------------------------------
nfolds = 5
nrepl = 10

# Data --------------------------------------------------------------------
ginv = G.inverse(G = Gmat, sparseform = TRUE)$Ginv.sparse
dat = transform(
  pheno,
  geno = as.factor(geno),
  env = as.factor(env),
  rn = as.factor(row_number),
  cn = as.factor(col_number),
  colgroup = as.factor(colgroup),
  rowgroup = as.factor(rowgroup),
  rc.group = as.factor(rc.group),
  pool = as.factor(pool),
  check = as.factor(check),
  genokeep = as.factor(ifelse(
    geno.cod %in% attr(ginv, 'rowNames'), geno.cod, NA
  )),
  genodrop = as.factor(ifelse(
    !geno.cod %in% attr(ginv, 'rowNames'), geno.cod, NA
  )),
  geno.cod = as.factor(geno.cod),
  check.cod = as.factor(check.cod)
)
str(dat)

load(file = "Saves/famod_sing.RDA")
modsel = lapply(famod_sing, expvar, env = 'env', geno = 'genokeep')

temp = data.frame(
  model = seq_along(famod_sing),
  ellapsed = do.call(c, lapply(famod_sing, function(x) x$ellapsed))
)

fasingsel = famod_sing[[2]]
rm(famod_sing)

name.env = levels(fasingsel$mf$env)
num.env = nlevels(fasingsel$mf$env)
name.geno = levels(fasingsel$mf$genokeep)
num.geno = nlevels(fasingsel$mf$genokeep)

vc = summary(fasingsel)$varcomp
fa = unique(na.exclude(str_extract(rownames(vc), "fa[0-9]")))
L = matrix(vc[grep('fa', rownames(vc)),1], nrow = num.env, ncol = length(fa),
           dimnames = list(name.env, fa))
Psi = diag(vc[which(grepl("vm", rownames(vc)) & !grepl("rr|var", rownames(vc))), 1])
dimnames(Psi) = list(name.env, name.env)
G = tcrossprod(L) + Psi
s2g_complete = diag(G)
initG = c(diag(Psi), L)

blu = summary(fasingsel, coef = TRUE)$coef.random
rr = data.frame(marg = blu[which(
  grepl("vm", rownames(blu)) &
    grepl("rr\\(env, i\\)", rownames(blu)) &
    !grepl("Comp", rownames(blu))
), 1])
rr = rr |> rownames_to_column('env') |>
  separate('env', into = c("env", 'geno'), sep = ':vm\\(genokeep, ginv\\)_') |>
  mutate(env = gsub("rr\\(env, i\\)_", '', env))
delta = data.frame(delta = blu[which(
  grepl("vm", rownames(blu)) &
    !grepl("rr\\(env, i\\)", rownames(blu)) &
    !grepl("Comp", rownames(blu))
), 1])
delta = delta |> rownames_to_column('env') |>
  separate('env', into = c("env", 'geno'), sep = ':vm\\(genokeep, ginv\\)_') |>
  mutate(env = gsub("env_", '', env))
ide = data.frame(nadd = blu[which(grepl("ide", rownames(blu))), 1]) |> 
  rownames_to_column('env') |>
  separate('env', into = c("env", 'geno'), sep = ':ide\\(genokeep\\)_') |>
  mutate(env = gsub("env_", '', env))

sinstgblup = left_join(rr, delta, by = c("geno", "env")) |> 
  left_join(ide, by = c("geno", 'env'))

blu = coef(fasingsel)$fixed
mu = blu[[1]]
pool = blu[grep('pool', rownames(blu))]; names(pool) = c("A", "B", "C")
env = blu[which(grepl("env", rownames(blu)) & !grepl("check", rownames(blu)))]

datcor = fasingsel$mf |> left_join(data.frame(
  pool = c("A", "B", "C"),
  pool_blue = pool
)) |> left_join(data.frame(
  env = gsub("env_", "", names(blu[which(grepl("env", rownames(blu)) & !grepl("check", rownames(blu))),])),
  env_blue = env
)) |> as.data.frame() |> 
  mutate(ycor = rytha - mu - pool_blue - env_blue)

### CV2 ---------------------------------------------------------------------
cvdata = list()
seed = 84
for (rept in 1:nrepl) {
  set.seed(seed * rept)
  cvdata[[rept]] = dat
  cvdata[[rept]]$set = NA
  for (id in levels(aux$geno.cod)) {
    cvdata[[rept]][which(cvdata[[rept]]$genokeep == id), 'set'] = sample(
      1:nfolds,
      size = dim(cvdata[[rept]][which(cvdata[[rept]]$genokeep == id), ])[1],
      replace = dim(cvdata[[rept]][which(cvdata[[rept]]$genokeep == id), ])[1] > nfolds
    )
  }
  cvdata[[rept]]$rept = rept
}

message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "initG", 'Gmat', "dat", "cvdata", "nfolds", "nrepl", 'datcor', 'sinstgblup', "s2g_complete"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("SS - CV2 // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = tryCatch({
      asreml(
        fixed = yNA  ~ pool + env*check.cod,
        random = ~ fa(env, 2, initG):vm(genokeep, ginv) + corh(env):ide(genokeep) +
          at(env):colgroup + at(env):rowgroup + at(env):colgroup:rowgroup,
        residual = ~ dsum(~ar1v(rn):ar1(cn)|env) ,
        sparse = ~ env:genodrop,
        na.action = na.method(x = "include"),
        data = w,
        maxit = 100,
        workspace = '5gb'
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      vc = summary(mod)$varcomp
      fa = unique(na.exclude(str_extract(rownames(vc), "fa[0-9]")))
      L = matrix(vc[grep('fa[0-9]', rownames(vc)),1], nrow = nlevels(w$env), ncol = length(fa),
                 dimnames = list(unique(w$env), fa))
      Psi = diag(vc[which(grepl("vm", rownames(vc)) & grepl("var", rownames(vc))), 1])
      dimnames(Psi) = list(unique(w$env), unique(w$env))
      G = tcrossprod(L) + Psi
      s2g_cv = diag(G)
      
      blu = summary(mod, coef = TRUE)$coef.random
      blup = as.data.frame(blu[which(
        grepl("vm", rownames(blu)) &
          # grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1:2])
      blup = blup |> rownames_to_column('env') |>
        separate('env', into = c("env", 'genokeep'), sep = ':vm\\(genokeep, ginv\\)_') |>
        mutate(env = gsub("fa\\(env, [0-9], initG\\)_", '', env)) 
      blup$reli_cv_F = NA  
      for (i in unique(blup$env)) {
        for (k in unique(blup$genokeep)) {
          blup[which(blup$env == i & blup$genokeep == k),'reli_cv_F'] = 1 - (blup[which(blup$env == i & blup$genokeep == k),'std.error']^2/(s2g_cv[which(names(s2g_cv)==i)]* diag(Gmat)[which(names(diag(Gmat)) == k)]))
        }
      }
      
      predvalu[[j]] = blup |> dplyr::select(-std.error) |> 
        left_join(unique(w[,c("geno", 'env', 'set', 'rept', 'rytha')]), 
                  by = c("genokeep" = 'geno', 'env')) |> 
        rename(y = rytha, ybar = solution, yhat = solution) |> 
        relocate('set', 'rept', .before = yhat)  |> 
        left_join(sinstgblup |> mutate(full_gebv = marg + delta) |> 
                    dplyr::select(geno, env, full_gebv), 
                  by = c('genokeep'='geno', 'env')) |> 
        filter(set == j) |> 
        left_join(unique(datcor[,c("genokeep", 'env', 'ycor')]), by = c('genokeep', 'env'))  
      
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/SS_CV2.RDS")

### CV1 ---------------------------------------------------------------------
set.seed(9684)
sets = split(
  rep(
    1:nfolds, length(unique(aux$geno.cod)) * nrepl
  )[order(runif(length(unique(aux$geno.cod)) * nrepl))],
  f = 1:nrepl
)
cvdata = lapply(sets, function(x){
  cvdata = dat
  cvdata = full_join(cvdata, data.frame(
    genokeep = unique(aux$geno.cod),
    set = x
  ), by = 'genokeep')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i

message("Starting cv")
cl = makeCluster(nrepl)
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics)})
clusterExport(cl, varlist = c("ginv", "initG", 'Gmat', "dat", "cvdata", "nfolds", "nrepl", 'datcor', 'sinstgblup', "s2g_complete"))
predval = parLapply(cl = cl, X = cvdata, fun = function(w){
  asreml.options(ai.sing = TRUE)
  predvalu = list()
  for (j in 1:nfolds) {
    message("SS - CV1 // fold: ", j)
    
    w$yNA = w$rytha
    w[which(w$set == j), c('yNA')] = NA
    
    mod = tryCatch({
      asreml(
        fixed = yNA  ~ pool + env*check.cod,
        random = ~ fa(env, 2, initG):vm(genokeep, ginv) + corh(env):ide(genokeep) +
          at(env):colgroup + at(env):rowgroup + at(env):colgroup:rowgroup,
        residual = ~ dsum(~ar1v(rn):ar1(cn)|env) ,
        sparse = ~ env:genodrop,
        na.action = na.method(x = "include"),
        data = w,
        maxit = 100,
        workspace = '5gb'
      )
    }, error = function(e){cat("The estimation was aborted or the model did not converge", 
                               fill = TRUE)})
    if(class(mod) == 'asreml' && mod$converge){
      vc = summary(mod)$varcomp
      fa = unique(na.exclude(str_extract(rownames(vc), "fa[0-9]")))
      L = matrix(vc[grep('fa[0-9]', rownames(vc)),1], nrow = nlevels(w$env), ncol = length(fa),
                 dimnames = list(unique(w$env), fa))
      Psi = diag(vc[which(grepl("vm", rownames(vc)) & grepl("var", rownames(vc))), 1])
      dimnames(Psi) = list(unique(w$env), unique(w$env))
      G = tcrossprod(L) + Psi
      s2g_cv = diag(G)
      
      blu = summary(mod, coef = TRUE)$coef.random
      blup = as.data.frame(blu[which(
        grepl("vm", rownames(blu)) &
          # grepl("rr\\(env, [0-9]\\)", rownames(blu)) &
          !grepl("Comp", rownames(blu))
      ), 1:2])
      blup = blup |> rownames_to_column('env') |>
        separate('env', into = c("env", 'genokeep'), sep = ':vm\\(genokeep, ginv\\)_') |>
        mutate(env = gsub("fa\\(env, [0-9], initG\\)_", '', env)) 
      blup$reli_cv = NA  
      for (i in unique(blup$env)) {
        for (k in unique(blup$genokeep)) {
          blup[which(blup$env == i & blup$genokeep == k),'reli_cv_F'] = 1 - (blup[which(blup$env == i & blup$genokeep == k),'std.error']^2/(s2g_cv[which(names(s2g_cv)==i)]* diag(Gmat)[which(names(diag(Gmat)) == k)]))
        }
      }

      predvalu[[j]] = blup |> dplyr::select(-std.error) |> 
        left_join(unique(w[,c("geno", 'env', 'set', 'rept', 'rytha')]), 
                  by = c("genokeep" = 'geno', 'env')) |> 
        rename(y = rytha, ybar = solution, yhat = solution) |> 
        relocate('set', 'rept', .before = yhat)  |> 
        left_join(sinstgblup |> mutate(full_gebv = marg + delta) |> 
                    dplyr::select(geno, env, full_gebv), 
                  by = c('genokeep'='geno', 'env')) |> 
        filter(set == j) |> 
        left_join(unique(datcor[,c("genokeep", 'env', 'ycor')]), by = c('genokeep', 'env'))  
      
    }else predvalu[[j]] = FALSE
  }
  predvalu = predvalu[which(do.call(c, lapply(predvalu, function(x)  !isFALSE(x))))]
  predvalu = do.call(rbind, predvalu)
  return(predvalu)
})

stopCluster(cl)

saveRDS(predval, file = "saves/resp2rev/SS_CV1_batch3.RDS")
