source("Codes/1_preparations.R")
library(lmmtools)

# First-stage analysis ----------------------------------------------------
data.list = split(pheno, pheno$env)

cl = makeCluster(length(data.list))
clusterEvalQ(cl, {library(asreml); library(tidyverse); library(ASRgenomics); library(lmmtools)})
clusterExport(cl, varlist = c("dereg", "up.mod", "Amat","pedigree", "Gmat", "data.list", "dat"))

stg1 = parLapply(cl = cl, X = data.list, fun = function(x){
  x = transform(
    x,
    geno = as.factor(geno),
    env = as.factor(env),
    rn = as.factor(row_number),
    cn = as.factor(col_number),
    colgroup = as.factor(colgroup),
    rowgroup = as.factor(rowgroup),
    rc.group = as.factor(rc.group),
    pool = as.factor(pool),
    check = as.factor(check),
    geno.cod = as.factor(geno.cod),
    check.cod = as.factor(check.cod)
  )
  x <<- x
  a = Sys.time()
  mod_fixed = stage1(
    fixed = rytha ~ pool + check.cod + geno.cod,
    random = ~ colgroup + rowgroup + colgroup:rowgroup,
    residual = ~ ar1v(rn):ar1(cn),
    na.action = na.method(x = "include"),
    data = x,
    type = "fixed",
    Genetic = "geno.cod",
    Trait = "env",
    n.trait = 1,
    s1trace = TRUE,
    pworkspace = '1gb',
    predGenetic = list(id = "geno.cod", classify = 'geno.cod:pool', rm = NULL),
    predict.options = list(ignore = "check.cod", levels = list(pool = c("A", "B")))
  )
  b = Sys.time()
  ella1 = b-a
  mod_fixed$pred.df$env = levels(x$env)
  
  ## Genotype as random (no kinship)
  mod_GI = asreml(
    fixed = rytha ~ pool + check.cod,
    random = ~ geno.cod + colgroup + rowgroup + colgroup:rowgroup,
    residual = ~ ar1v(rn):ar1(cn),
    na.action = na.method(x = "include"),
    data = x,
    maxit = 50,
    workspace = '1gb'
  )
  if (any(na.exclude(mod_GI$vparameters.pc > 1))) mod_GI = up.mod(mod_GI)

  
  vc.GI = as.data.frame(summary(mod_GI)$varcomp) |>
    rownames_to_column("effect") |>
    mutate(
      model = "BLUP",
      trait = 'rytha',
      effect = ifelse(grepl("geno", effect), 'geno', effect)
    )
  h2.GI = vpredict(mod_GI, H2 ~ V4 / (V4 + V7)) |>
    mutate(herita = "H2",
           entry = "BLUP",
           trait = 'rytha')
  
  aux = na.omit(unique(x[,c("geno.cod",'pool')]))
  
  pred = predict(
    mod_GI,
    classify = "geno.cod:pool",
    vcov = TRUE,
    pworkspace = "2gb", ignore = "check.cod", levels = list(pool = c("A", "B"))
  )
  pred$vcov = pred$vcov[which(
    paste(pred$pvals$geno.cod, pred$pvals$pool, sep = '@') %in%
      paste(aux$geno.cod, aux$pool, sep = '@')
  ), which(
    paste(pred$pvals$geno.cod, pred$pvals$pool, sep = '@') %in%
      paste(aux$geno.cod, aux$pool, sep = '@')
  )]
  pred$pvals = pred$pvals[which(
    paste(pred$pvals$geno.cod, pred$pvals$pool, sep = '@') %in%
      paste(aux$geno.cod, aux$pool, sep = '@')
  ), ]
  
  if (any(pred$pvals$status == "Aliased")) {
    pred$vcov = pred$vcov[-which(pred$pvals$status == "Aliased"), -which(pred$pvals$status == "Aliased")]
    pred$pvals = pred$pvals[-which(pred$pvals$status == "Aliased"), ]
  }
  blu = summary(mod_GI, coef = TRUE)$coef.random
  blup = pred$pvals |>
    mutate(
      env = unique(x$env),
      reli = 1 - (std.error^2 / vc.GI[grepl("geno", vc.GI$effect), 'component'])
    ) |>
    dplyr::select(-status, -std.error) |>
    # filter(grepl("UGP", geno.cod)) |> 
    rename(geno = geno.cod, mublup = predicted.value) |>
    left_join(
      data.frame(blu[which(grepl("geno", rownames(blu))), ]) |>
        rownames_to_column("geno") |>
        mutate_at("geno", str_replace, ".*geno\\.cod\\_", "") |>
        dplyr::select(-z.ratio, -std.error) |>
        rename(blup = solution), by = 'geno'
    )
  
  lrtest.GI = as.data.frame(lrt(
    mod_GI,
    asreml(
      fixed = rytha  ~ pool + check.cod,
      random = ~ colgroup + rowgroup + colgroup:rowgroup,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    )
  ), boundary = TRUE) |>
    mutate(entry = "BLUP", trait = "rytha")
  a = Sys.time()
  mod_GI = stage1(
    fixed = rytha ~ pool + check.cod,
    random = ~ geno.cod + colgroup + rowgroup + colgroup:rowgroup,
    residual = ~ ar1v(rn):ar1(cn),
    na.action = na.method(x = "include"),
    data = x,
    type = "random",
    Genetic = "geno.cod",
    Trait = "env",
    n.trait = 1,
    s1trace = TRUE,
    pworkspace = '1gb',
    predGenetic = list(id = "geno.cod", classify = 'geno.cod:pool', rm = NULL),
    predict.options = list(ignore = "check.cod", levels = list(pool = c("A", "B")))
  )
  b = Sys.time()
  ella2 = b-a
  mod_GI$pred.df$env = levels(x$env)
  
  ### Genotype as random (animal model)
  Asparse = ASRgenomics::full2sparse(Amat)
  attr(Asparse, "INVERSE") = FALSE
  Asparse <<- Asparse
  x$genod = as.factor(x$geno.cod)
  x <<- x
  a = Sys.time()
  mod_GAA = stage1(
    fixed = rytha ~ pool + check.cod,
    random = ~ vm(geno.cod, Asparse) + ide(geno.cod) + colgroup + rowgroup + colgroup:rowgroup,
    residual = ~ ar1v(rn):ar1(cn),
    na.action = na.method(x = "include"),
    data = x,
    type = "random",
    Genetic = "geno.cod",
    Trait = "env",
    n.trait = 1,
    s1trace = TRUE,
    pworkspace = '1gb',
    predGenetic = list(id = "geno.cod", classify = 'geno.cod:pool', rm = NULL),
    predict.options = list(ignore = "check.cod", levels = list(pool = c("A", "B")))
  )
  b = Sys.time()
  ella2 = b-a
  mod_GA$pred.df$env = levels(x$env)
  
  
  mod_GA = asreml(
    fixed = rytha  ~ pool + check.cod,
    random = ~ vm(geno.cod, Asparse) + ide(geno.cod) + colgroup + rowgroup + colgroup:rowgroup,
    residual = ~ ar1v(rn):ar1(cn),
    na.action = na.method(x = "include", y = "include"),
    data = x,
    maxit = 50,
    workspace = '1gb'
  )
  if (any(na.exclude(mod_GA$vparameters.pc > 1)))
    mod_GA = up.mod(mod_GA)

  Kbar = mean(Amat)
  diagK = mean(diag(Amat))
  Dk <<- diagK-Kbar
  
  vc.GA = as.data.frame(summary(mod_GA)$varcomp) |>
    rownames_to_column("effect") |>
    mutate(
      model = "ABLUP",
      D = Dk,
      trait = "rytha",
      effect = case_when(
        grepl("vm", effect) ~ "add",
        grepl("ide", effect) ~ "nadd",
        .default = effect
      )
    )

  h2.GA = rbind(
    vpredict(mod_GA, H2 ~ Dk*V4 / (V4 + V5 + V8)) |>
      mutate(
        herita = "h2",
        entry = "ABLUP",
        trait = "rytha"
      ),
    vpredict(mod_GA, H2 ~ (Dk*V4 + V5) / (V4 + V5 + V8)) |>
      mutate(
        herita = "H2",
        entry = "ABLUP",
        trait = 'rytha'
      )
  )
  
  add = as.data.frame(lrt.asreml(
    asreml(
      fixed = rytha  ~ pool + check.cod,
      random = ~ vm(geno.cod, Asparse) + colgroup + rowgroup + colgroup:rowgroup,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    ),
    asreml(
      fixed = rytha  ~ pool + check.cod,
      random = ~ geno.cod + colgroup + rowgroup + colgroup:rowgroup,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    )
  ), boundary = TRUE) |> mutate(entry = "ABLUP", trait = 'rytha')
  nadd = as.data.frame(lrt(
    asreml(
      fixed = rytha  ~ pool + check.cod,
      random = ~ vm(geno.cod, Asparse) + ide(geno.cod) +
        colgroup + rowgroup + colgroup:rowgroup,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    ),
    asreml(
      fixed = rytha  ~ pool + check.cod,
      random = ~ vm(geno.cod, Asparse) + colgroup + rowgroup + colgroup:rowgroup,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    )
  ), boundary = TRUE) |> mutate(entry = "ABLUP", trait = 'rytha')
  lrtest.GA = list(add = add, nadd = nadd)
  
  pred = predict(
    mod_GA,
    classify = "geno.cod:pool",
    vcov = TRUE,
    pworkspace = "2gb", ignore = "check.cod", levels = list(pool = c("A", "B"))
  )
  pred$vcov = pred$vcov[which(
    paste(pred$pvals$geno.cod, pred$pvals$pool, sep = '@') %in%
      paste(aux$geno.cod, aux$pool, sep = '@')
  ), which(
    paste(pred$pvals$geno.cod, pred$pvals$pool, sep = '@') %in%
      paste(aux$geno.cod, aux$pool, sep = '@')
  )]
  pred$pvals = pred$pvals[which(
    paste(pred$pvals$geno.cod, pred$pvals$pool, sep = '@') %in%
      paste(aux$geno.cod, aux$pool, sep = '@')
  ), ]
  
  if (any(pred$pvals$status == "Aliased")) {
    pred$vcov = pred$vcov[-which(pred$pvals$status == "Aliased"), -which(pred$pvals$status == "Aliased")]
    pred$pvals = pred$pvals[-which(pred$pvals$status == "Aliased"), ]
  }
  
  blu = summary(mod_GA, coef = TRUE)$coef.random
  ablup = data.frame(blu[which(grepl("vm", rownames(blu))), ]) |>
    rownames_to_column("geno") |>
    mutate_at("geno", str_replace, ".*\\)\\_", "") |>
    dplyr::select(-z.ratio) |>
    rename(ablup = solution) |>
    mutate(
      ablup = ifelse(ablup == 0, NA, ablup),
      std.error = ifelse(is.na(ablup), NA, std.error),
      reli = 1 - (std.error^2 / vc.GA[grepl("^add", vc.GA$effect), 'component']),
      env = unique(x$env)
    ) |> left_join(pred$pvals[,c("geno.cod", 'predicted.value')] |> 
                     rename(tgv = predicted.value), by = c("geno"= 'geno.cod'))
  
  list(
    models = list(
      BLUE = mod_fixed,
      BLUP = mod_GI,
      ABLUP = mod_GAA
    ),
    selection = list(
      BLUP = blup,
      ABLUP = ablup
    ),
    vc = list(BLUP = vc.GI, ABLUP = vc.GA),
    herita = list(h2.GI, h2.GA),
    # cvalsel = predat,
    lrtest = list(lrtest.GI, lrtest.GA),
    ellapsed = list(
      BLUE = ella1, 
      BLUP = ella2, 
      ABLUP = ella3
    )
  )
  
})
stopCluster(cl)

models = lapply(stg1, function(x) x$models)
selection = lapply(stg1, function(x) x$selection)
vc = lapply(stg1, function(x) x$vc)
herita = lapply(stg1, function(x) x$herita)
lrtest =lapply(stg1, function(x) x$lrtest)
ellapsed = lapply(stg1, function(x) x$ellapsed)

save(models, file = "Saves/stg1_models.RDA")
save(selection, file = "Saves/stg1_selection.RDA")
save(vc, file = "Saves/stg1_vc.RDA")
save(herita, file = "Saves/stg1_herita.RDA")
save(lrtest, file = "Saves/stg1_lrtest.RDA")
save(ellapsed, file = "Saves/stg1_ellapsed.RDA")
