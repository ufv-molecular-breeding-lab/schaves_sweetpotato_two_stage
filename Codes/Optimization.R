rm(list=ls())

# Loading the preparations ------------------------------------------------
source("Codes/Preparations.R")
library(parallel)
library(TSDFGS)

# Preparing the marker data -----------------------------------------------
snp = snp[which(rownames(snp) %in% rownames(Gmat)),]
snp = snp[-duplicated(rownames(snp)),]
snpA = snp[which(rownames(snp) %in% unique(pheno[which(pheno$pool == "A"), c('geno')])),]
snpB = snp[which(rownames(snp) %in% unique(pheno[which(pheno$pool == "B"), c('geno')])),]

snp_scale = scale(snp, center = TRUE, scale = FALSE)
snpA_scale = scale(snpA, center = TRUE, scale = FALSE)
snpB_scale = scale(snpB, center = TRUE, scale = FALSE)

## PCA --------------
MlM = svd(snpA_scale)
L_A = snpA_scale %*% MlM$v
MlM = svd(snpB_scale)
L_B = snpB_scale %*% MlM$v

# Optimization ------------------------------------------------------------
trsize_A = ceiling(nrow(L_A) * seq(.5, .8, .05))
trsize_B = ceiling(nrow(L_B) * seq(.5, .8, .05))

utg.optA = mclapply(
  mc.cores = 41,
  X = 1:40,
  FUN = function(x) {
    optsize = list()
    for (i in 1:length(trsize_A)) {
      rscore = optTrain(
        L_A,
        cand = 1:nrow(L_A),
        n.train = trsize_A[i],
        method = 'rScore',
        console = FALSE
      )
      pev = optTrain(
        L_A,
        cand = 1:nrow(L_A),
        n.train = trsize_A[i],
        method = 'PEV',
        console = FALSE
      )
      cd = optTrain(
        L_A,
        cand = 1:nrow(L_A),
        n.train = trsize_A[i],
        method = 'CD',
        console = FALSE
      )
      optsize[[i]] = list(rscore = rscore,
                          pev = pev,
                          cd = cd)
    }
    names(optsize) = paste0('size_', trsize_A)
    optsize
  }
)
save(utg.optA, file = "Saves/trt_opt_poolA.RDA")


utg.optB = mclapply(
  mc.cores = 41,
  X = 1:40,
  FUN = function(x) {
    optsize = list()
    for (i in 1:length(trsize_B)) {
      rscore = optTrain(
        L_B,
        cand = 1:nrow(L_B),
        n.train = trsize_B[i],
        method = 'rScore',
        console = FALSE
      )
      pev = optTrain(
        L_B,
        cand = 1:nrow(L_B),
        n.train = trsize_B[i],
        method = 'PEV',
        console = FALSE
      )
      cd = optTrain(
        L_B,
        cand = 1:nrow(L_B),
        n.train = trsize_B[i],
        method = 'CD',
        console = FALSE
      )
      optsize[[i]] = list(rscore = rscore,
                          pev = pev,
                          cd = cd)
    }
    names(optsize) = paste0('size_', trsize_B)
    optsize
  }
)
save(utg.optB, file = "Saves/trt_opt_poolB.RDA")

# Processing --------------------------------------------------------------

## Pool A ----------------
names(optA) = paste0("R", seq_along(optA))
optdfA = do.call(rbind,
                 lapply(optA, function(x)
                   do.call(
                     rbind,
                     lapply(x, function(x)
                       do.call(rbind, lapply(x, function(x)
                         data.frame(geno = rownames(
                           L_A
                         )[x$OPTtrain]))) |> rownames_to_column("crit") |>
                         mutate_at("crit", str_replace, '\\..*', ''))
                   ) |> rownames_to_column("trt_size") |>
                     mutate_at("trt_size", str_replace, '\\..*', ''))) |> rownames_to_column('repl') |>
  mutate_at('repl', str_replace, '\\..*', '')

freqlistA = lapply(split(optdfA, optdfA$trt_size), function(x) as.data.frame(table(x$geno, x$crit)))
freqlistA = lapply(freqlistA, function(x)
  lapply(split(x, x$Var2), function(x)
    x |> mutate(Freq = ifelse(Var1 %in% rownames(L_A)[duplicated(rownames(L_A))], 
                              Freq/2, Freq)) |>  arrange(desc(Freq))))
for(i in names(freqlistA)){
  for(j in names(freqlistA[[1]])){
    aux = rep(0, nrow(freqlistA[[i]][[j]]))
    aux[1:as.numeric(gsub('size_','',i))] = 1
    freqlistA[[i]][[j]][,'trtopt'] = aux
  }
}

freqdf = data.frame(geno = rownames(L_A)) |>
  full_join(
    do.call(rbind, lapply(freqlistA, function(x) do.call(rbind, x))) |> 
      rownames_to_column('trt_size') |>
      mutate_at("trt_size", str_replace, '\\..*', ''),
    by = c('geno' = "Var1")
  )

pca.A = FactoMineR::PCA(Gmat[rownames(Gmat) %in% levels(entry_A$geno),
                             colnames(Gmat) %in% levels(entry_A$geno)], 
                        scale.unit = TRUE, graph = FALSE)
pc.df.A = data.frame(
  geno = names(pca.A$ind$coord[,1]),
  pc1 = pca.A$ind$coord[,1], 
  pc2 = pca.A$ind$coord[,2], row.names = NULL
) |> full_join(freqdf, by = 'geno') |> filter(!is.na(Freq)) |> 
  mutate(smooth = Freq/max(Freq))

saveRDS(pc.df.A, file = 'Saves/optA.rds')

## Pool B ----------------
load("saves/optimization_step/trt_opt_poolB.RDA")

lapply(
  lapply(optB, function(x)
    do.call(rbind, lapply(x, function(x)
      data.frame(diff = do.call(rbind, lapply(x, function(x) {
        x$TOPscore[length(x$TOPscore)] - x$TOPscore[length(x$TOPscore) - 1]
      }))) |> rownames_to_column("crit"))) |> rownames_to_column("trt_size") |> 
      mutate_at("trt_size", str_replace, '\\..*','')),
  function(x) all(x$diff == 0)
)

names(optB) = paste0("R", seq_along(optB))
optdfB = do.call(rbind,
                 lapply(optB, function(x)
                   do.call(
                     rbind,
                     lapply(x, function(x)
                       do.call(rbind, lapply(x, function(x)
                         data.frame(geno = rownames(
                           L_B
                         )[x$OPTtrain]))) |> rownames_to_column("crit") |>
                         mutate_at("crit", str_replace, '\\..*', ''))
                   ) |> rownames_to_column("trt_size") |>
                     mutate_at("trt_size", str_replace, '\\..*', ''))) |> rownames_to_column('repl') |>
  mutate_at('repl', str_replace, '\\..*', '')

freqlistB = lapply(split(optdfB, optdfB$trt_size), function(x) as.data.frame(table(x$geno, x$crit)))
freqlistB = lapply(freqlistB, function(x)
  lapply(split(x, x$Var2), function(x)
    x |> mutate(Freq = ifelse(Var1 %in% rownames(L_B)[duplicated(rownames(L_B))], 
                              Freq/2, Freq)) |>  arrange(desc(Freq))))
for(i in names(freqlistB)){
  for(j in names(freqlistB[[1]])){
    aux = rep(0, nrow(freqlistB[[i]][[j]]))
    aux[1:as.numeric(gsub('size_','',i))] = 1
    freqlistB[[i]][[j]][,'trtopt'] = aux
  }
}

freqdf = data.frame(geno = rownames(L_B)) |>
  full_join(
    do.call(rbind, lapply(freqlistB, function(x) do.call(rbind, x))) |> 
      rownames_to_column('trt_size') |>
      mutate_at("trt_size", str_replace, '\\..*', ''),
    by = c('geno' = "Var1")
  )

pca.B = FactoMineR::PCA(Gmat[rownames(Gmat) %in% levels(entry_B$geno),
                             colnames(Gmat) %in% levels(entry_B$geno)], 
                        scale.unit = TRUE, graph = FALSE)
pc.df.B = data.frame(
  geno = names(pca.B$ind$coord[,1]),
  pc1 = pca.B$ind$coord[,1], 
  pc2 = pca.B$ind$coord[,2], row.names = NULL
) |> full_join(freqdf, by = 'geno') |> filter(!is.na(Freq)) |> 
  mutate(smooth = Freq/max(Freq))

saveRDS(pc.df.B, file = 'Saves/optB.rds')

# Cross-validation --------------------------------------------------------

## Pool A ------------------------------------------------------------------
pc.df.A = readRDS(file = 'Saves/optA.rds')
optset = lapply(split(pc.df.A, pc.df.A$trt_size), function(x){
  lapply(split(x, x$Var2), function(x) unique(x[which(x$trtopt == 1), 'geno']))
})

opt_cv = lapply(optset, function(x) lapply(x,function(x){
  asreml.options(ai.sing = TRUE, maxit = 50, workspace = '2gb', pworkspace = '2gb')
  
  df = ablup[which(ablup$geno %in% data[which(data$pool == "A"),"geno"]),]
  df$yNA = df[, "dAblup_adj"]
  df$yNA[!df$geno %in% x] = NA
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  
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
  
  prdablup_adj = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
    rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
    separate("geno", into = c("geno", "env"), sep = ":env_") %>%
    mutate_at("geno", str_replace, ".*_", "") %>%
    mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
    left_join(df[,c("geno","env","dAblup_adj",'yNA')]) %>% 
    left_join(blue[,c('geno','env','blue')]) |> 
    filter(!geno %in% x) %>% rename(gebv = solution, ybar = blue)
  
  prdablup_adj
  
}))
save(opt_cv,file = 'Saves/opt_cv_A.RDA')

opt_cv_resA = do.call(rbind, lapply(opt_cv, function(x) do.call(rbind, lapply(x, function(x){
  data.frame(
    corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
    mspe = mean((x$ybar - x$gebv)^2, na.rm = TRUE)
  )
})) %>% rownames_to_column("crit"))) %>% rownames_to_column('trt_size') %>% 
  mutate_at('trt_size',str_replace, '\\..*','') |> 
  mutate_at("trt_size", str_replace, 'size_','')


sizes = as.numeric(gsub('size_','',unique(pc.df.A$trt_size)))
random_cv = mclapply(X = sizes, mc.cores = length(sizes), function(x){
  asreml.options(ai.sing = TRUE, maxit = 50, workspace = '2gb', pworkspace = '2gb')
  i = 1
  df = ablup[which(ablup$geno %in% data[which(data$pool == "A"),"geno"]),]
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  result = list()
  repeat
  {
    set.seed(8 * i)
    ran_train = sample(levels(df$geno), x)
    df$yNA = df[, "dAblup_adj"]
    df$yNA[!df$geno %in% ran_train] = NA
    
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
    
    prdablup_adj = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","dAblup_adj", 'yNA')]) %>% 
      left_join(blue[,c('geno','env','blue')])  %>% 
      filter(!geno %in% ran_train) %>% rename(gebv = solution, ybar = blue)
    
    result[[i]] = data.frame(
      corr = cor(prdablup_adj$gebv, prdablup_adj$ybar, use = 'na.or.complete'),
      mspe = mean((prdablup_adj$ybar - prdablup_adj$gebv)^2, na.rm = TRUE)
    )
    
    i = i + 1
    if (i > 100)
      break
  }
  do.call(rbind, result) %>% mutate(size = x)
})

save(random_cv, file = 'Saves/random_cv_A.RDA')


## Pool B ------------------------------------------------------------------
pc.df.B = readRDS(file = 'Saves/optB.rds')
optset = lapply(split(pc.df.B, pc.df.B$trt_size), function(x){
  lapply(split(x, x$Var2), function(x) unique(x[which(x$trtopt == 1), 'geno']))
})

opt_cv = lapply(optset, function(x) lapply(x,function(x){
  asreml.options(ai.sing = TRUE, maxit = 50, workspace = '2gb', pworkspace = '2gb')
  
  df = ablup[which(ablup$geno %in% data[which(data$pool == "B"),"geno"]),]
  df$yNA = df[, "dAblup_adj"]
  df$yNA[!df$geno %in% x] = NA
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  
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
  
  prdablup_adj = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
    rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
    separate("geno", into = c("geno", "env"), sep = ":env_") %>%
    mutate_at("geno", str_replace, ".*_", "") %>%
    mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
    left_join(df[,c("geno","env","dAblup_adj")]) %>% 
    left_join(blue[,c('geno','env','blue')]) |> 
    filter(!geno %in% x) %>% rename(gebv = solution, ybar = blue)
  
  prdablup_adj
  
}))
save(opt_cv, file = 'Saves/opt_cv_B.RDA')

opt_cv_resB = do.call(rbind, lapply(opt_cv, function(x) do.call(rbind, lapply(x, function(x){
  data.frame(
    corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
    mspe = mean((x$ybar - x$gebv)^2, na.rm = TRUE)
  )
})) %>% rownames_to_column("crit"))) %>% rownames_to_column('trt_size') %>% 
  mutate_at('trt_size',str_replace, '\\..*','') |> 
  mutate_at("trt_size", str_replace, 'size_','')


sizes = as.numeric(gsub('size_','',unique(pc.df.B$trt_size)))
random_cv = mclapply(X = sizes, mc.cores = length(sizes), function(x){
  asreml.options(ai.sing = TRUE, maxit = 50, workspace = '2gb', pworkspace = '2gb')
  i = 1
  df = ablup[which(ablup$geno %in% data[which(data$pool == "B"),"geno"]),]
  df = transform(df, geno = as.factor(geno), env = as.factor(env))
  df = droplevels(df[-which(is.na(df$WABLUP_adj)), ])
  temp = Gmat[rownames(Gmat) %in% levels(df$geno), 
              colnames(Gmat) %in% levels(df$geno)]
  Gblup = full2sparse(temp)
  attr(Gblup, 'INVERSE') = FALSE
  result = list()
  repeat
  {
    set.seed(8 * i)
    ran_train = sample(levels(df$geno), x)
    df$yNA = df[, "dAblup_adj"]
    df$yNA[!df$geno %in% ran_train] = NA
    
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
    
    prdablup_adj = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
      rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
      separate("geno", into = c("geno", "env"), sep = ":env_") %>%
      mutate_at("geno", str_replace, ".*_", "") %>%
      mutate(entry = "dABLUP_adj", trait = 'rytha', converge = mod$converge) %>% 
      left_join(df[,c("geno","env","dAblup_adj", 'yNA')]) %>% 
      left_join(blue[,c('geno','env','blue')])  %>% 
      filter(!geno %in% ran_train) %>% rename(gebv = solution, ybar = blue)
    
    result[[i]] = data.frame(
      corr = cor(prdablup_adj$gebv, prdablup_adj$ybar, use = 'na.or.complete'),
      mspe = mean((prdablup_adj$ybar - prdablup_adj$gebv)^2, na.rm = TRUE)
    )
    
    i = i + 1
    if (i > 100)
      break
  }
  do.call(rbind, result) %>% mutate(size = x)
})

save(random_cv, file = 'Saves/random_cv_B.RDA')



