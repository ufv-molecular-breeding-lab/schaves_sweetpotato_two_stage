rm(list=ls())

# Loading the preparations ------------------------------------------------
source("Codes/Preparations.R")
load(file = "Saves/RES_1st.RDA")
library(parallel)
library(ASRgenomics)

# Preparations for each scenario ------------------------------------------
ablup = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind, x$entries[which(names(x$entries) == 'ABLUP')])
  rownames(aux) = NULL
  aux
})) |> rownames_to_column('env') |> mutate_at("env", str_replace, "\\..*",'') |> 
  mutate(entry = "ABLUP") |> rename(WABLUP = WBLUP) |> 
  filter(geno %in% rownames(Gmat)) |> droplevels() |> 
  left_join(unique(pheno[,c("geno",'pool')]))
ablup[which(ablup$reli < 0.15),'WABLUP'] = NA
ablup[which(ablup$reli_adj < 0.15),c("dAblup_adj",'WABLUP_adj')] = NA

blup = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind, x$entries[which(names(x$entries) == 'BLUP')])
  rownames(aux) = NULL
  aux
})) |> rownames_to_column('env') |> mutate_at("env", str_replace, "\\..*",'') |> 
  filter(geno %in% rownames(Gmat)) |> droplevels() |> 
  left_join(unique(pheno[,c("geno",'pool')]))

blue = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind, x$entries[which(names(x$entries) == 'BLUE')])
  rownames(aux) = NULL
  droplevels(aux)
})) |> rownames_to_column('env') |> mutate_at("env", str_replace, "\\..*",'') |> 
  filter(geno %in% rownames(Gmat)) |> droplevels()

# Cross-validation II: within pool assessment ---------------------------------
nfolds = 5 
nrept = 5
seed = 78

cvdata2 = list(
  BLUE = list(),
  dBLUP = list(),
  dABLUP = list(),
  dABLUP_adj = list()
)

for (k in names(cvdata2)) {
  if(k == "BLUE") temp = droplevels(blue[which(blue$pool == "A"),])
  else if (k == "dBLUP") temp = droplevels(blup[which(blup$pool == "A"),])
  else if (k == "dABLUP" | k == 'dABLUP_adj') temp = droplevels(ablup[which(ablup$pool == "A"),])
  for (j in 1:nrept) {
    set.seed(seed+j)
    cvdata2[[k]][[j]] = temp
    cvdata2[[k]][[j]]$set = NA
    for (i in unique(temp$geno)) {
      cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i, 'set'] = sample(
        1:nfolds, 
        size = dim(cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i,])[1],
        replace = dim(cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i,])[1] > nfolds
      )
    }
    cvdata2[[k]][[j]]['repl'] = j
  }
}

cv2.poolA.blue = mclapply(X = cvdata2$BLUE, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblue = full2sparse(temp)
  attr(Gblue, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("BLUE - CV2 (A) // fold:", fold)
    
    df$yNA = df[, "blue"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblue):corh(env) + ide(geno),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUE",
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "BLUE", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blue")]) %>%
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv2.poolA.blup = mclapply(X = cvdata2$dBLUP, mc.cores = nrept, function(df){
  preds = list()
  if(any(is.na(df$WBLUP))) df = df[-which(is.na(df$WBLUP)), ]
  df = transform(df, geno = as.factor(geno), env = as.factor(env), genod = as.factor(geno))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dBLUP - CV2 (A) // fold:", fold)
    
    df$yNA = df[, "dblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env) + idv(genod),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUP",
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dBLUP", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blup")]) %>%
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv2.poolA.dAblup = mclapply(X = cvdata2$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP - CV2 (A) // fold:", fold)
    df$yNA = df[, "dAblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP",
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv2.poolA.dAblup_adj = mclapply(X = cvdata2$dABLUP_adj, mc.cores = nrept, function(df){
  preds = list()
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP_adj - CV2 (A) // fold:", fold)
    
    df$yNA = df[, "dAblup_adj"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP_adj",
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup_adj")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv2.poolA = list(
  BLUE = cv2.poolA.blue,
  dBLUP = cv2.poolA.blup,
  dABLUP = cv2.poolA.dAblup,
  dABLUP_adj = cv2.poolA.dAblup_adj
)

cvdata2 = list(
  BLUE = list(),
  dBLUP = list(),
  dABLUP = list(),
  dABLUP_adj = list()
)

for (k in names(cvdata2)) {
  if(k == "BLUE") temp = droplevels(blue[which(blue$pool == "B"),])
  else if (k == "dBLUP") temp = droplevels(blup[which(blup$pool == "B"),])
  else if (k == "dABLUP" | k == 'dABLUP_adj') temp = droplevels(ablup[which(ablup$pool == "B"),])
  for (j in 1:nrept) {
    set.seed(seed+j)
    cvdata2[[k]][[j]] = temp
    cvdata2[[k]][[j]]$set = NA
    for (i in unique(temp$geno)) {
      cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i, 'set'] = sample(
        1:nfolds, 
        size = dim(cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i,])[1],
        replace = dim(cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i,])[1] > nfolds
      )
    }
    cvdata2[[k]][[j]]['repl'] = j
  }
}


cv2.poolB.blue = mclapply(X = cvdata2$BLUE, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblue = full2sparse(temp)
  attr(Gblue, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("BLUE - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "blue"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblue):corh(env) + ide(geno),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUE",
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "BLUE", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blue")]) %>%
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv2.poolB.blup = mclapply(X = cvdata2$dBLUP, mc.cores = nrept, function(df){
  preds = list()
  df = df[-which(is.na(df$WBLUP)), ]
  df = transform(df, geno = as.factor(geno), env = as.factor(env), genod = as.factor(geno))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dBLUP - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "dblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env) + idv(genod),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUP",
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dBLUP", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blup")]) %>%
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv2.poolB.dAblup = mclapply(X = cvdata2$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "dAblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP",
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv2.poolB.dAblup_adj = mclapply(X = cvdata2$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP_adj - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "dAblup_adj"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP_adj",
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup_adj")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv2.poolB = list(
  BLUE = cv2.poolB.blue,
  dBLUP = cv2.poolB.blup,
  dABLUP = cv2.poolB.dAblup,
  dABLUP_adj = cv2.poolB.dAblup_adj
)

cv2 = list(A = cv2.poolA, B = cv2.poolB)

save(cv2, file = "Saves/cv2_intrapool.RDA")

# Cross-validation II: using both pools -----------------------------------
nfolds = 5 
nrept = 5
seed = 1997

cvdata2 = list(
  BLUE = list(),
  dBLUP = list(),
  dABLUP = list(),
  dABLUP_adj = list()
)

for (k in names(cvdata2)) {
  if(k == "BLUE") temp = blue
  else if (k == "dBLUP") temp = blup
  else if (k == "dABLUP" | k == 'dABLUP_adj') temp = ablup
  for (j in 1:nrept) {
    set.seed(seed+j)
    cvdata2[[k]][[j]] = temp
    cvdata2[[k]][[j]]$set = NA
    for (i in unique(temp$geno)) {
      cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i, 'set'] = sample(
        1:nfolds, 
        size = dim(cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i,])[1],
        replace = dim(cvdata2[[k]][[j]][cvdata2[[k]][[j]]$geno == i,])[1] > nfolds
      )
    }
    cvdata2[[k]][[j]]['repl'] = j
  }
}

cv2.blue = mclapply(X = cvdata2$BLUE, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblue = full2sparse(temp)
  attr(Gblue, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("BLUE - CV2 (B) // fold:", fold)
    
    
    df$yNA = df[, "blue"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblue):corh(env) + ide(geno),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUE",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "BLUE", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blue")]) %>%
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv2.blup = mclapply(X = cvdata2$dBLUP, mc.cores = nrept, function(df){
  preds = list()
  
  df = df[-which(is.na(df$WBLUP)), ]
  df = transform(df, geno = as.factor(geno), env = as.factor(env), genod = as.factor(geno))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dBLUP - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "dblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env) + idv(genod),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUP",
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dBLUP", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blup")]) %>%
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv2.dAblup = mclapply(X = cvdata2$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "dAblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv2.dAblup_adj = mclapply(X = cvdata2$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP_adj - CV2 (B) // fold:", fold)
    
    df$yNA = df[, "dAblup_adj"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP_adj",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup_adj")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv2 = list(
  BLUE = cv2.blue,
  dBLUP = cv2.blup,
  dABLUP = cv2.dAblup,
  dABLUP_adj = cv2.dAblup_adj
)

save(cv2, file = "Saves/cv2_twopool.RDA")

# Cross-validation I: within pool assessment ---------------------------------------------------------
nfolds = 5 
nrept = 5
seed = 5156

cvdata = list(
  BLUE = list(),
  dBLUP = list(),
  dABLUP = list(),
  dABLUP_adj = list()
)

for (k in names(cvdata)) {
  if(k == "BLUE") temp = droplevels(blue[which(blue$pool == "A"),])
  else if (k == "dBLUP") temp = droplevels(blup[which(blup$pool == "A"),])
  else if (k == "dABLUP" | k == 'dABLUP_adj') temp = droplevels(ablup[which(ablup$pool == "A"),])
  set.seed(seed)
  sets = split(
    rep(
      1:nfolds, length(unique(temp$geno)) * nrept
    )[order(runif(length(unique(temp$geno)) * nrept))],
    f = 1:nrept
  )
  cvdata[[k]] = lapply(sets, function(x){
    cvdata = temp
    cvdata = merge(cvdata, data.frame(
      geno = unique(temp$geno),
      set = x
    ), by = 'geno')
  })
  
  for (i in 1:length(cvdata[[k]])) cvdata[[k]][[i]]$repl = i
}

cv1.poolA.blue = mclapply(X = cvdata$BLUE, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblue = full2sparse(temp)
  attr(Gblue, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("BLUE - CV1 (A) // fold:", fold)
    
    df$yNA = df[, "blue"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblue):corh(env) + ide(geno),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUE",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "BLUE", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blue")]) %>%
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv1.poolA.blup = mclapply(X = cvdata$dBLUP, mc.cores = nrept, function(df){
  preds = list()
  df = df[-which(is.na(df$WBLUP)), ]
  df = transform(df, geno = as.factor(geno), env = as.factor(env), genod = as.factor(geno))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dBLUP - CV1 (A) // fold:", fold)
    
    df$yNA = df[, "dblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env) + idv(genod),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUP",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dBLUP", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blup")]) %>%
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv1.poolA.dAblup = mclapply(X = cvdata$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP - CV1 (A) // fold:", fold)
    
    df$yNA = df[, "dAblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv1.poolA.dAblup_adj = mclapply(X = cvdata$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP_adj - CV1 (A) // fold:", fold)
    
    df$yNA = df[, "dAblup_adj"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP_adj",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup_adj")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv1.poolA = list(
  BLUE = cv1.poolA.blue,
  dBLUP = cv1.poolA.blup,
  dABLUP = cv1.poolA.dAblup,
  dABLUP_adj = cv1.poolA.dAblup_adj
)

cvdata = list(
  BLUE = list(),
  dBLUP = list(),
  dABLUP = list(),
  dABLUP_adj = list()
)

for (k in names(cvdata)) {
  if(k == "BLUE") temp = droplevels(blue[which(blue$pool == "B"),])
  else if (k == "dBLUP") temp = droplevels(blup[which(blup$pool == "B"),])
  else if (k == "dABLUP" | k == 'dABLUP_adj') temp = droplevels(ablup[which(ablup$pool == "B"),])
  set.seed(seed)
  sets = split(
    rep(
      1:nfolds, length(unique(temp$geno)) * nrept
    )[order(runif(length(unique(temp$geno)) * nrept))],
    f = 1:nrept
  )
  cvdata[[k]] = lapply(sets, function(x){
    cvdata = temp
    cvdata = merge(cvdata, data.frame(
      geno = unique(temp$geno),
      set = x
    ), by = 'geno')
  })
  
  for (i in 1:length(cvdata[[k]])) cvdata[[k]][[i]]$repl = i
}

cv1.poolB.blue = mclapply(X = cvdata$BLUE, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblue = full2sparse(temp)
  attr(Gblue, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("BLUE - CV1 (B) // fold:", fold)
    
    df$yNA = df[, "blue"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblue):corh(env) + ide(geno),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUE",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "BLUE", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blue")]) %>%
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv1.poolB.blup = mclapply(X = cvdata$dBLUP, mc.cores = nrept, function(df){
  preds = list()
  df = df[-which(is.na(df$WBLUP)), ]
  df = transform(df, geno = as.factor(geno), env = as.factor(env), genod = as.factor(geno))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dBLUP - CV1 (B) // fold:", fold)
    
    
    df$yNA = df[, "dblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env) + idv(genod),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUP",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dBLUP", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blup")]) %>%
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv1.poolB.dAblup = mclapply(X = cvdata$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP - CV1 (B) // fold:", fold)
    
    df$yNA = df[, "dAblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv1.poolB.dAblup_adj = mclapply(X = cvdata$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP_adj - CV1 (B) // fold:", fold)
    
    df$yNA = df[, "dAblup_adj"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP_adj",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup_adj")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv1.poolB = list(
  BLUE = cv1.poolB.blue,
  dBLUP = cv1.poolB.blup,
  dABLUP = cv1.poolB.dAblup,
  dABLUP_adj = cv1.poolB.dAblup_adj
)

cv1 = list(A = cv1.poolA, B = cv1.poolB)

save(cv1, file = "Saves/cv1_intrapool.RDA")

# Cross-validation I: using both pools ---------------------------------------------------------
nfolds = 5 
nrept = 5
seed = 985

for (k in names(cvdata)) {
  if(k == "BLUE") temp = blue
  else if (k == "dBLUP") temp = blup
  else if (k == "dABLUP" | k == 'dABLUP_adj') temp = ablup
  set.seed(seed)
  sets = split(
    rep(
      1:nfolds, length(unique(temp$geno)) * nrept
    )[order(runif(length(unique(temp$geno)) * nrept))],
    f = 1:nrept
  )
  cvdata[[k]] = lapply(sets, function(x){
    cvdata = temp
    cvdata = merge(cvdata, data.frame(
      geno = unique(temp$geno),
      set = x
    ), by = 'geno')
  })
  
  for (i in 1:length(cvdata[[k]])) cvdata[[k]][[i]]$repl = i
}

cv1.blue = mclapply(X = cvdata$BLUE, mc.cores = nrept, function(df){
  preds = list()
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblue = full2sparse(temp)
  attr(Gblue, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("BLUE - CV1 // fold:", fold)
    
    
    df$yNA = df[, "blue"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblue):corh(env) + ide(geno),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUE",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "BLUE", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blue")]) %>%
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv1.blup = mclapply(X = cvdata$dBLUP, mc.cores = nrept, function(df){
  preds = list()
  if(any(is.na(df$WBLUP))) df = df[-which(is.na(df$WBLUP)), ]
  df = transform(df, geno = as.factor(geno), env = as.factor(env), genod = as.factor(geno))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dBLUP - CV1 // fold:", fold)
    
    df$yNA = df[, "dblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env) + idv(genod),
      family = asr_gaussian(dispersion = 1),
      weights = "WBLUP",
      data = df,
      workspace = '1.5gb',
      maxit = 50,
      ai.sing = TRUE
    )
    
    blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
    preds[[fold]] = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dBLUP", trait = 'rytha', converge = mod$converge) |> 
      left_join(df[,c("geno","env","set","blup")]) %>%
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
  }
  do.call(rbind, preds)
})

cv1.dAblup = mclapply(X = cvdata$dABLUP, mc.cores = nrept, function(df){
  preds = list()
  if(any(is.na(df$WABLUP)))  df = droplevels(df[-which(is.na(df$WABLUP)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP - CV1 // fold:", fold)
    
    df$yNA = df[, "dAblup"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP",
      # na.action = na.method(x = "include", y = "include"),
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv1.dAblup_adj = mclapply(X = cvdata$dABLUP_adj, mc.cores = nrept, function(df){
  preds = list()
  if(any(is.na(df$WABLUP_adj)))  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  for (fold in 1:nfolds) {
    message("dABLUP_adj - CV1 // fold:", fold)
    
    df$yNA = df[, "dAblup_adj"]
    df$yNA[df$set == fold] = NA
    
    mod = asreml(
      fixed = yNA ~ env,
      random = ~ vm(geno, Gblup):corh(env),
      family = asr_gaussian(dispersion = 1),
      weights = "WABLUP_adj",
      data = df,
      workspace = '1.5gb',
      maxit = 50
    )
    
    preds[[fold]] = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","set","dAblup_adj")]) %>% 
      left_join(blue[,c('geno','env','blue')]) |> 
      filter(set == fold) %>% rename(gebv = solution, ybar = blue)
    
  }
  do.call(rbind, preds)
})

cv1 = list(
  BLUE = cv1.blue,
  dBLUP = cv1.blup,
  dABLUP = cv1.dAblup,
  dABLUP_adj = cv1.dAblup_adj
)

save(cv1, file = "Saves/cv1_twopool.RDA")