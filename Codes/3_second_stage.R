source("Codes/1_preparations.R")
load(file = "Saves/stg1_models.RDA")

# BLUE -------------------------------------------------------------------

## Two pools (diagonal and full weights) ---------------------------------
blue = do.call(rbind, lapply(models, function(x) x$BLUE$pred.df))
V_blue = do.call(bdiag, lapply(models, function(x) x$BLUE$bdVt))
blue = droplevels(blue[which(grepl("UGP", blue$geno.cod)), ])
ginv = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blue$geno.cod),
                          colnames(Gmat) %in% levels(blue$geno.cod)],
                 sparseform = TRUE)$Ginv.sparse
V_blue = V_blue[which(blue$geno.cod %in% attr(ginv, 'rowNames')), 
                which(blue$geno.cod %in% attr(ginv, 'rowNames'))]
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

famodblue_dw = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_diag),
    data = blue,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodblue_dw[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_diag),
    data = blue,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 30
  )
  if(any(na.exclude(famodblue_dw[[paste0("FA", i)]]$vparameters.pc > 1))) famodblue_dw[[paste0("FA", i)]] = up.mod(famodblue_dw[[paste0("FA", i)]])
  b = Sys.time()
  famodblue_dw[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblue_dw[[paste0("FA", i)]]$converge) break
}
save(famodblue_dw, file = "saves/resp2rev/famodblue_dw.RDA")

famodblue_fw = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
    data = blue,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodblue_fw[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
    data = blue,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 30
  )
  if(any(na.exclude(famodblue_fw[[paste0("FA", i)]]$vparameters.pc > 1))) famodblue_fw[[paste0("FA", i)]] = up.mod(famodblue_fw[[paste0("FA", i)]])
  b = Sys.time()
  famodblue_fw[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblue_fw[[paste0("FA", i)]]$converge) break
}
save(famodblue_fw, file = "saves/resp2rev/famodblue_fw.RDA")

## Per pool (diagonal and full weights) ---------------------------------
V_blueA = V_blue[which(blue$pool == "A"), which(blue$pool == "A")]
blueA = droplevels(blue[which(blue$pool == "A"), ])
ginvA = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blueA$geno.cod), colnames(Gmat) %in% levels(blueA$geno.cod)], sparseform = TRUE)$Ginv.sparse

V_diagA = matrix(0, nrow = nrow(V_blueA), ncol = ncol(V_blueA))
diag(V_diagA) = 1/diag(solve(V_blueA))
rownames(V_diagA) = rownames(V_blueA)
colnames(V_diagA) = colnames(V_blueA)
attr(V_diagA, 'INVERSE') = FALSE

V_blueB = V_blue[which(blue$pool == "B"), which(blue$pool == "B")]
blueB = droplevels(blue[which(blue$pool == "B"), ])
ginvB = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blueB$geno.cod), colnames(Gmat) %in% levels(blueB$geno.cod)], sparseform = TRUE)$Ginv.sparse

V_diagB = matrix(0, nrow = nrow(V_blueB), ncol = ncol(V_blueB))
diag(V_diagB) = 1/diag(solve(V_blueB))
rownames(V_diagB) = rownames(V_blueB)
colnames(V_diagB) = colnames(V_blueB)
attr(V_diagB, 'INVERSE') = FALSE

famodblue_dw_A = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagA),
    data = blueA,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodblue_dw_A[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagA),
    data = blueA,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 50
  )
  if(any(na.exclude(famodblue_dw_A[[paste0("FA", i)]]$vparameters.pc > 1))) famodblue_dw_A[[paste0("FA", i)]] = up.mod(famodblue_dw_A[[paste0("FA", i)]])
  b = Sys.time()
  famodblue_dw_A[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblue_dw_A[[paste0("FA", i)]]$converge) break
}

famodblue_dw_B = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagB),
    data = blueB,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodblue_dw_B[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagB),
    data = blueB,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 50
  )
  if(any(na.exclude(famodblue_dw_B[[paste0("FA", i)]]$vparameters.pc > 1))) famodblue_dw_B[[paste0("FA", i)]] = up.mod(famodblue_dw_B[[paste0("FA", i)]])
  b = Sys.time()
  famodblue_dw_B[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblue_dw_B[[paste0("FA", i)]]$converge) break
}
save(famodblue_dw_A, famodblue_dw_B, file = "saves/resp2rev/famodblue_dw_pool.RDA")


V_blue_sparse = full2sparse(as.matrix(V_blueA))
attr(V_blue_sparse, 'INVERSE') = FALSE
famodblue_fw_A = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
    data = blueA,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodblue_fw_A[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
    data = blueA,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 50
  )
  if(any(na.exclude(famodblue_fw_A[[paste0("FA", i)]]$vparameters.pc > 1))) famodblue_fw_A[[paste0("FA", i)]] = up.mod(famodblue_fw_A[[paste0("FA", i)]])
  b = Sys.time()
  famodblue_fw_A[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblue_fw_A[[paste0("FA", i)]]$converge) break
}

V_blue_sparse = full2sparse(as.matrix(V_blueB))
attr(V_blue_sparse, 'INVERSE') = FALSE
famodblue_fw_B = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
    data = blueB,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodblue_fw_B[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_blue_sparse),
    data = blueB,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 50
  )
  if(any(na.exclude(famodblue_fw_B[[paste0("FA", i)]]$vparameters.pc > 1))) famodblue_fw_B[[paste0("FA", i)]] = up.mod(famodblue_fw_B[[paste0("FA", i)]])
  b = Sys.time()
  famodblue_fw_B[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblue_fw_B[[paste0("FA", i)]]$converge) break
}

save(famodblue_fw_A,famodblue_fw_B, file = "saves/resp2rev/famodblue_fw_pool.RDA")

# dBLUP -------------------------------------------------------------------
## Two pools (diagonal and full weights) ---------------------------------
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

famodblup_dw = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_diag),
    data = blup,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodblup_dw[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_diag),
    data = blup,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 30
  )
  
  if(any(na.exclude(famodblup_dw[[paste0("FA", i)]]$vparameters.pc > 1))) famodblup_dw[[paste0("FA", i)]] = up.mod(famodblup_dw[[paste0("FA", i)]])
  b = Sys.time()
  famodblup_dw[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblup_dw[[paste0("FA", i)]]$converge) break
}
save(famodblup_dw, file = "saves/resp2rev/famodblup_dw.RDA")

famodblup_fw = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
    data = blup,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodblup_fw[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
    data = blup,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 30
  )
  if(any(na.exclude(famodblup_fw[[paste0("FA", i)]]$vparameters.pc > 1))) famodblup_fw[[paste0("FA", i)]] = up.mod(famodblup_fw[[paste0("FA", i)]])
  b = Sys.time()
  famodblup_fw[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblup_fw[[paste0("FA", i)]]$converge) break
}
save(famodblup_fw, file = "saves/resp2rev/famodblup_fw.RDA")

## Per pool (diagonal and full weights) ---------------------------------
V_blupA = V_blup[which(blup$pool == "A"), which(blup$pool == "A")]
blupA = droplevels(blup[which(blup$pool == "A"), ])
ginvA = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blupA$geno.cod), colnames(Gmat) %in% levels(blupA$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_blupB = V_blup[which(blup$pool == "B"), which(blup$pool == "B")]
blupB = droplevels(blup[which(blup$pool == "B"), ])
ginvB = G.inverse(G = Gmat[rownames(Gmat) %in% levels(blupB$geno.cod), colnames(Gmat) %in% levels(blupB$geno.cod)], sparseform = TRUE)$Ginv.sparse


V_diagA = matrix(0, nrow = nrow(V_blupA), ncol = ncol(V_blupA))
diag(V_diagA) = 1/diag(solve(V_blupA))
rownames(V_diagA) = rownames(V_blupA)
colnames(V_diagA) = colnames(V_blupA)
attr(V_diagA, 'INVERSE') = FALSE

famodblup_dw_A = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagA),
    data = blupA,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodblup_dw_A[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagA),
    data = blupA,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 50
  )
  if(any(na.exclude(famodblup_dw_A[[paste0("FA", i)]]$vparameters.pc > 1))) famodblup_dw_A[[paste0("FA", i)]] = up.mod(famodblup_dw_A[[paste0("FA", i)]])
  b = Sys.time()
  famodblup_dw_A[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblup_dw_A[[paste0("FA", i)]]$converge) break
}

V_diagB = matrix(0, nrow = nrow(V_blupB), ncol = ncol(V_blupB))
diag(V_diagB) = 1/diag(solve(V_blupB))
rownames(V_diagB) = rownames(V_blupB)
colnames(V_diagB) = colnames(V_blupB)
attr(V_diagB, 'INVERSE') = FALSE
famodblup_dw_B = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagB),
    data = blupB,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodblup_dw_B[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagB),
    data = blupB,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 50
  )
  
  if(any(na.exclude(famodblup_dw_B[[paste0("FA", i)]]$vparameters.pc > 1))) famodblup_dw_B[[paste0("FA", i)]] = up.mod(famodblup_dw_B[[paste0("FA", i)]])
  b = Sys.time()
  famodblup_dw_B[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblup_dw_B[[paste0("FA", i)]]$converge) break
}
save(famodblup_dw_A, famodblup_dw_B, file = "saves/resp2rev/famodblup_dw_pool.RDA")


V_blup_sparse = full2sparse(as.matrix(V_blupA))
attr(V_blup_sparse, 'INVERSE') = FALSE
famodblup_fw_A = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
    data = blupA,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodblup_fw_A[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
    data = blupA,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 30
  )
  if(any(na.exclude(famodblup_fw_A[[paste0("FA", i)]]$vparameters.pc > 1))) famodblup_fw_A[[paste0("FA", i)]] = up.mod(famodblup_fw_A[[paste0("FA", i)]])
  b = Sys.time()
  famodblup_fw_A[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblup_fw_A[[paste0("FA", i)]]$converge) break
}

V_blup_sparse = full2sparse(as.matrix(V_blupB))
attr(V_blup_sparse, 'INVERSE') = FALSE
famodblup_fw_B = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
    data = blupB,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodblup_fw_B[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_blup_sparse),
    data = blupB,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 30
  )
  if(any(na.exclude(famodblup_fw_B[[paste0("FA", i)]]$vparameters.pc > 1))) famodblup_fw_B[[paste0("FA", i)]] = up.mod(famodblup_fw_B[[paste0("FA", i)]])
  b = Sys.time()
  famodblup_fw_B[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodblup_fw_B[[paste0("FA", i)]]$converge) break
}

save(famodblup_fw_A,famodblup_fw_B, file = "saves/resp2rev/famodblup_fw_pool.RDA")


# dABLUP -------------------------------------------------------------------
## Two pools (diagonal and full weights) ---------------------------------
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

famodablup_dw = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_diag),
    data = ablup,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodablup_dw[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv) + 
      corh(env):ide(geno.cod) + vm(TG, V_diag),
    data = ablup,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 30
  )
  
  if(any(na.exclude(famodablup_dw[[paste0("FA", i)]]$vparameters.pc > 1))) famodablup_dw[[paste0("FA", i)]] = up.mod(famodablup_dw[[paste0("FA", i)]])
  b = Sys.time()
  famodablup_dw[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodablup_dw[[paste0("FA", i)]]$converge) break
}
save(famodablup_dw, file = "saves/resp2rev/famodablup_dw.RDA")

famodablup_fw = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv)+
      corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
    data = ablup,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodablup_fw[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginv) + diag(env):vm(geno.cod, ginv)+
      corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
    data = ablup,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 30
  )
  if(any(na.exclude(famodablup_fw[[paste0("FA", i)]]$vparameters.pc > 1))) famodablup_fw[[paste0("FA", i)]] = up.mod(famodablup_fw[[paste0("FA", i)]])
  b = Sys.time()
  famodablup_fw[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodablup_fw[[paste0("FA", i)]]$converge) break
}
save(famodablup_fw, file = "saves/resp2rev/famodablup_fw.RDA")

## Per pool (diagonal and full weights) ---------------------------------
V_ablupA = V_ablup[which(ablup$pool == "A"), which(ablup$pool == "A")]
ablupA = droplevels(ablup[which(ablup$pool == "A"), ])
ginvA = G.inverse(G = Gmat[rownames(Gmat) %in% levels(ablupA$geno.cod), colnames(Gmat) %in% levels(ablupA$geno.cod)], sparseform = TRUE)$Ginv.sparse
V_ablupB = V_ablup[which(ablup$pool == "B"), which(ablup$pool == "B")]
ablupB = droplevels(ablup[which(ablup$pool == "B"), ])
ginvB = G.inverse(G = Gmat[rownames(Gmat) %in% levels(ablupB$geno.cod), colnames(Gmat) %in% levels(ablupB$geno.cod)], sparseform = TRUE)$Ginv.sparse

V_diagA = matrix(0, nrow = nrow(V_ablupA), ncol = ncol(V_ablupA))
diag(V_diagA) = 1/diag(solve(V_ablupA))
rownames(V_diagA) = rownames(V_ablupA)
colnames(V_diagA) = colnames(V_ablupA)
attr(V_diagA, 'INVERSE') = FALSE

famodablup_dw_A = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagA),
    data = ablupA,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodablup_dw_A[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagA),
    data = ablupA,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 50
  )
  
  if(any(na.exclude(famodablup_dw_A[[paste0("FA", i)]]$vparameters.pc > 1))) famodablup_dw_A[[paste0("FA", i)]] = up.mod(famodablup_dw_A[[paste0("FA", i)]])
  b = Sys.time()
  famodablup_dw_A[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodablup_dw_A[[paste0("FA", i)]]$converge) break
}

V_diagB = matrix(0, nrow = nrow(V_ablupB), ncol = ncol(V_ablupB))
diag(V_diagB) = 1/diag(solve(V_ablupB))
rownames(V_diagB) = rownames(V_ablupB)
colnames(V_diagB) = colnames(V_ablupB)
attr(V_diagB, 'INVERSE') = FALSE
famodablup_dw_B = list()
for(i in 1:4){
  a = Sys.time()
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagB),
    data = ablupB,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  famodablup_dw_B[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB) + 
      corh(env):ide(geno.cod) + vm(TG, V_diagB),
    data = ablupB,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "3gb",
    maxit = 50
  )
  if(any(na.exclude(famodablup_dw_B[[paste0("FA", i)]]$vparameters.pc > 1))) famodablup_dw_B[[paste0("FA", i)]] = up.mod(famodablup_dw_B[[paste0("FA", i)]])
  b = Sys.time()
  famodablup_dw_B[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodablup_dw_B[[paste0("FA", i)]]$converge) break
}
save(famodablup_dw_A, famodablup_dw_B, file = "saves/resp2rev/famodablup_dw_pool.RDA")


V_ablup_sparse = full2sparse(as.matrix(V_ablupA))
attr(V_ablup_sparse, 'INVERSE') = FALSE
famodablup_fw_A = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA)+
      corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
    data = ablupA,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodablup_fw_A[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvA) + diag(env):vm(geno.cod, ginvA)+
      corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
    data = ablupA,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 50
  )
  if(any(na.exclude(famodablup_fw_A[[paste0("FA", i)]]$vparameters.pc > 1))) famodablup_fw_A[[paste0("FA", i)]] = up.mod(famodablup_fw_A[[paste0("FA", i)]])
  b = Sys.time()
  famodablup_fw_A[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodablup_fw_A[[paste0("FA", i)]]$converge) break
}

V_ablup_sparse = full2sparse(as.matrix(V_ablupB))
attr(V_ablup_sparse, 'INVERSE') = FALSE
famodablup_fw_B = list()
for(i in 1:4){
  mod = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB)+
      corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
    data = ablupB,
    family = asr_gaussian(dispersion = 0.0001), start.values = TRUE
  )
  params = mod$vparameters.table
  params$Value[grep("TG", params$Component)] = 1
  params$Constraint[grep("TG", params$Component)] = "F"
  
  a = Sys.time()
  famodablup_fw_B[[paste0("FA", i)]] = asreml(
    fixed = rytha ~ -1 + env,
    random = ~ rr(env, i):vm(geno.cod, ginvB) + diag(env):vm(geno.cod, ginvB)+
      corh(env):ide(geno.cod) + vm(TG, V_ablup_sparse),
    data = ablupB,
    family = asr_gaussian(dispersion = 0.0001),
    G.param = params,
    workspace = "5gb",
    maxit = 50
  )
  if(any(na.exclude(famodablup_fw_B[[paste0("FA", i)]]$vparameters.pc > 1))) famodablup_fw_B[[paste0("FA", i)]] = up.mod(famodablup_fw_B[[paste0("FA", i)]])
  b = Sys.time()
  famodablup_fw_B[[paste0("FA", i)]]$ellapsed = b - a
  if(!famodablup_fw_B[[paste0("FA", i)]]$converge) break
}

save(famodablup_fw_A,famodablup_fw_B, file = "saves/resp2rev/famodablup_fw_pool.RDA")
