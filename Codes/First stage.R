rm(list=ls())

# Loading the preparations ------------------------------------------------
source("Codes/Preparations.R")

# First-stage analysis ----------------------------------------------------
data.list = split(pheno, pheno$env)
data.list = data.list[which(grepl("OT", names(data.list)))]

## This process will run in parallel
results = parallel::mclapply(
  X = data.list,
  mc.cores = (length(data.list) + 1),  # Adjust the number of cores here
  FUN = function(x) {
  
    x = transform(
      x,
      geno = as.factor(geno),
      rn = as.factor(row_number),
      cn = as.factor(col_number),
      colgroup = as.factor(colgroup),
      rowgroup = as.factor(rowgroup),
      rc.group = as.factor(rc.group),
      pool = as.factor(ifelse(pool == "C", NA, pool)),
      check = as.factor(ifelse(status == "Check", "Check", "Candidate")),
      geno.cod = as.factor(ifelse(status != "Check", geno, NA)),
      check.cod = as.factor(ifelse(status == "Check", geno, NA))
    )
    
    x <<- x
    
    ## Genotype as fixed
    modf = asreml(
      fixed = rytha ~ -1 + pool +
        check.cod+
        geno.cod,
      random = ~ rn + cn + colgroup + rowgroup + rc.group,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    )
    if (any(na.exclude(modf$vparameters.pc > 1)))
      modf = gencomp:::.up.mod(modf)
    
    pred = predict(
      modf,
      classify = "geno.cod:pool",
      vcov = TRUE,
      pworkspace = "2gb", ignore = "check.cod"
    )
    
    
    if (any(pred$pvals$status == "Aliased")) {
      pred$vcov = pred$vcov[-which(pred$pvals$status == "Aliased"), 
                            -which(pred$pvals$status == "Aliased")]
      pred$pvals = pred$pvals[-which(pred$pvals$status == "Aliased"), ]
    }
    
    blue = pred$pvals |> select(-status, -std.error) |>
      rename(blue = predicted.value, geno = geno.cod) |>
      mutate(WBLUE = diag(solve(pred$vcov)),
             entry = "BLUE",
             trait = "rytha") |>
      filter(grepl("UGP", geno))
    
    ## Genotype as random (no kinship)
    mod_GI = asreml(
      fixed = rytha  ~ pool +
        check.cod,
      random = ~ geno.cod +
        rn + cn + colgroup +
        rowgroup + rc.group,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    )
    if (any(na.exclude(mod_GI$vparameters.pc > 1)))
      mod_GI = up.mod(mod_GI)
    
    vc.GI = as.data.frame(summary(mod_GI)$varcomp) |>
      rownames_to_column("effect") |>
      mutate(
        model = "BLUP",
        trait = 'rytha',
        effect = ifelse(grepl("geno", effect), 'geno', effect)
      )
    h2.GI = vpredict(mod_GI, H2 ~ V6 / (V6 + V9)) |>
      mutate(herita = "H2",
             entry = "BLUP",
             trait = 'rytha')
    
    
    pred = predict(mod_GI,
                   classify = "geno.cod",
                   vcov = TRUE,
                   pworkspace = "2gb", ignore = "check.cod")
    
    if (any(pred$pvals$status == "Aliased")) {
      pred$vcov = pred$vcov[-which(pred$pvals$status == "Aliased" ), 
                            -which(pred$pvals$status == "Aliased")]
      pred$pvals = pred$pvals[-which(pred$pvals$status == "Aliased"), ]
    }
    
    blu = summary(mod_GI, coef = TRUE)$coef.random
    blup = data.frame(blu[which(grepl("geno", rownames(blu))), ]) |>
      rownames_to_column("geno") |>
      mutate_at("geno", str_replace, ".*geno\\.cod\\_", "") |>
      select(-z.ratio) |>
      filter(grepl("UGP", geno)) |> rename(blup = solution) |>
      mutate(
        blup = ifelse(blup == 0, NA, blup),
        std.error = ifelse(is.na(blup), NA, std.error),
        reli = 1 - std.error ^ 2 / vc.GI[grepl("geno", vc.GI$effect), 'component'],
        dblup = ifelse(reli < 0.15, NA, blup / reli)
      ) |>
      left_join(
        pred$pvals |>
          rename(mublup = predicted.value, geno = geno.cod) |>
          mutate(WBLUP = diag(solve(pred$vcov)), entry = "BLUP") |>
          select(-std.error, -status) |>
          filter(grepl("UGP", geno)),
        by = 'geno'
      ) |>
      mutate(
        mublup = ifelse(is.na(blup), NA, mublup),
        WBLUP = ifelse(is.na(blup), NA, WBLUP),
        dmublup = ifelse(reli < 0.15, NA, mublup /
                           reli),
        trait = 'rytha'
      )
    
    lrtest.GI = as.data.frame(lrt(
      mod_GI,
      asreml(
        fixed = rytha  ~ pool+
          check.cod,
        random = ~ rn + cn + colgroup +
          rowgroup + rc.group,
        residual = ~ ar1v(rn):ar1(cn),
        na.action = na.method(x = "include", y = "include"),
        data = x,
        maxit = 50,
        workspace = '1gb'
      )
    )) |>
      mutate(entry = "BLUP", trait = "rytha")
    
    ## Genotype as random (Amatrix)
    Asparse = ASRgenomics::full2sparse(Amat[rownames(Amat) %in% unique(c(as.character(x$geno), x$male, x$female)), colnames(Amat) %in% unique(c(as.character(x$geno), x$male, x$female))])
    attr(Asparse, "INVERSE") = FALSE
    
    Asparse <<- Asparse
    
    x$genod = x$geno.cod
    
    mod_GA = asreml(
      fixed = rytha  ~ pool+
        check.cod,
      random = ~ vm(geno.cod, Asparse) +
        genod+
        rn + cn + colgroup +
        rowgroup + rc.group,
      residual = ~ ar1v(rn):ar1(cn),
      na.action = na.method(x = "include", y = "include"),
      data = x,
      maxit = 50,
      workspace = '1gb'
    )
    if (any(na.exclude(mod_GA$vparameters.pc > 1)))
      mod_GA = up.mod(mod_GA)
    
    
    vc.GA = as.data.frame(summary(mod_GA)$varcomp) |>
      rownames_to_column("effect") |>
      mutate(
        model = "ABLUP",
        trait = "rytha",
        effect = case_when(
          grepl("vm", effect) ~ "genoadd",
          grepl("genod", effect) ~ "genod",
          .default = effect
        )
      )
    h2.GA = rbind(
      vpredict(mod_GA, H2 ~ V6 / (V6 + V7 + V10)) |>
        mutate(
          herita = "h2",
          entry = "ABLUP",
          trait = "rytha"
        ),
      vpredict(mod_GA, H2 ~ (V6 + V7) / (V6 + V7 + V10)) |>
        mutate(
          herita = "H2",
          entry = "ABLUP",
          trait = 'rytha'
        )
    )
    
    lrtest.GA = list(
      A = as.data.frame(lrt(
        mod_GA,
        asreml(
          fixed = rytha  ~ pool +
            check.cod,
          random = ~ genod +
            rn + cn + colgroup +
            rowgroup + rc.group,
          residual = ~ ar1v(rn):ar1(cn),
          na.action = na.method(x = "include", y = "include"),
          data = x,
          maxit = 50,
          workspace = '1gb'
        )
      )) |> mutate(entry = "ABLUP", trait = 'rytha'),
      `NA` = as.data.frame(lrt(
        mod_GA,
        asreml(
          fixed = rytha  ~ pool +
            check.cod,
          random = ~ vm(geno.cod, Asparse) +
            rn + cn + colgroup +
            rowgroup + rc.group,
          residual = ~ ar1v(rn):ar1(cn),
          na.action = na.method(x = "include", y = "include"),
          data = x,
          maxit = 50,
          workspace = '1gb'
        )
      )) |> mutate(entry = "ABLUP", trait = 'rytha')
    )
    lrtest.GA = data.frame(do.call(rbind, lrtest.GA)) |>
      rownames_to_column("effect")
    
    pred = predict(mod_GA,
                   classify = "geno.cod",
                   vcov = TRUE,
                   pworkspace = "2gb",
                   ignore = "check.cod")
    
    if (any(pred$pvals$status == "Aliased")) {
      pred$vcov = pred$vcov[-which(pred$pvals$status == "Aliased"), 
                            -which(pred$pvals$status == "Aliased" )]
      pred$pvals = pred$pvals[-which(pred$pvals$status == "Aliased"), ]
    }
    
    blu = summary(mod_GA, coef = TRUE)$coef.random
    ablup = data.frame(blu[which(grepl("vm", rownames(blu))), ]) |>
      rownames_to_column("geno") |>
      mutate_at("geno", str_replace, ".*\\)\\_", "") |>
      select(-z.ratio) |>
      rename(Ablup = solution) |>
      mutate(
        Ablup = ifelse(Ablup == 0, NA, Ablup),
        std.error = ifelse(is.na(Ablup), NA, std.error),
        reli = 1 - std.error ^ 2 / vc.GA[grepl("genoadd", vc.GA$effect), 'component'],
        dAblup = ifelse(reli < 0.15, NA, Ablup / reli)
      ) |>
      left_join(
        pred$pvals |>
          rename(muAblup = predicted.value, geno = geno.cod) |>
          mutate(WBLUP = diag(solve(pred$vcov)), entry = "ABLUP") |>
          select(-std.error, -status),
        by = 'geno'
      ) |>
      mutate(
        muAblup = ifelse(is.na(Ablup), NA, muAblup),
        WBLUP = ifelse(is.na(Ablup), NA, WBLUP),
        dmuAblup = ifelse(reli < 0.15, NA, muAblup /
                            reli),
        trait = 'rytha'
      )
    
    ## Garrick et al. (2009) deregression
    blupar = ablup[which(!grepl("UGP", ablup$geno)), ]
    ablup = ablup[which(grepl("UGP", ablup$geno)), ]
    
    lambda = c((1 - h2.GA[which(h2.GA$herita == "h2"), 'Estimate']) /
                 h2.GA[which(h2.GA$herita == "h2"), 'Estimate'])
    ablup[, c("dAblup_adj", "reli_adj")] = NA
    for (g in ablup$geno) {
      gi = ablup[which(ablup$geno == g), "Ablup"]
      ri = ablup[which(ablup$geno == g), "reli"]
      gpa = mean(blupar[which(blupar$geno %in%
                                ped[which(ped$geno == g), c("male", "female")]), "Ablup"])
      rpa = sum(blupar[which(blupar$geno %in%
                               ped[which(ped$geno == g), c("male", "female")]), "reli"]) /
        4
      alfa = 1 / (0.5 - rpa)
      cpapa = (0.5 - rpa) / lambda
      cii = (1 - ri) / lambda
      delta = cpapa / cii
      ZpaZpa = lambda * (0.5 * alfa - 4) + 0.5 * lambda * sqrt((alfa ^
                                                                  2 + 16) / delta)
      ZiZi = delta * ZpaZpa + 2 * lambda * (2 * delta - 1)
      
      rhs = matrix(
        c(ZpaZpa + 4 * lambda, -2 * lambda, -2 * lambda, ZiZi + 2 * lambda),
        nrow = 2,
        ncol = 2,
        byrow = TRUE
      ) %*% matrix(c(gpa, gi))
      gamma_star = rhs[2]
      dblup_adj = gamma_star / ZiZi
      ablup[which(ablup$geno == g), "dAblup_adj"] = dblup_adj
      reli_adj = 1 - lambda / (ZiZi + lambda)
      ablup[which(ablup$geno == g), "reli_adj"] = reli_adj
    }
    
    ## Cross-validation for defining the best c value 
    Gsparse = ASRgenomics::full2sparse(Gmat[rownames(Gmat) %in% ablup$geno, colnames(Gmat) %in% ablup$geno])
    attr(Gmat, 'INVERSE') = FALSE

    temp.blu = ablup[which(ablup$geno %in% rownames(Gmat)), ]
    temp.blu$geno = as.factor(temp.blu$geno)
    cv.cval = list()
    for (k in 1:length(seq(0.05, 0.95, 0.05))) {
      c.val = seq(0.05, 0.95, 0.05)[k]
      temp.blu = temp.blu |>
        mutate(w_adj = (1 - c(h2.GA[which(h2.GA$herita == "h2"), 'Estimate'])) /
                 (c.val + (1 - reli_adj) /
                    reli_adj) * c(h2.GA[which(h2.GA$herita == "h2"), 'Estimate']))

      nrepl = 5
      nfolds = 5
      seed = 13

      cvdata = list()
      for (rept in 1:nrepl) {
        set.seed(seed * rept)
        cvdata[[rept]] = temp.blu[, c('geno', 'dAblup_adj', 'w_adj')]
        cvdata[[rept]]$set = NA
        for (id in temp.blu$geno) {
          cvdata[[rept]][cvdata[[rept]]$geno == id, 'set'] = sample(
            1:nfolds,
            size = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1],
            replace = dim(cvdata[[rept]][cvdata[[rept]]$geno == id, ])[1] > nfolds
          )
        }
      }

      for (u in 1:length(cvdata))
        cvdata[[u]]$repl = u

      cv = lapply(cvdata, function(df2) {
        res.list = list()
        cvstats = list()
        for (fold in 1:nfolds) {
          df2$yNA = df2$dAblup_adj
          df2$yNA[df2$set == fold] = NA

          mod = asreml(
            fixed =  yNA ~ 1,
            random = ~ vm(geno, Gsparse),
            weights = 'w_adj',
            family = asr_gaussian(dispersion = 1),
            na.action = na.method(x = "include", y = "include"),
            data = df2,
            maxit = 50,
            workspace = "1gb"
          )
          pred = predict(mod, classify = 'geno')$pvals
          pred = left_join(pred[, 1:2], df2[, c("geno", 'set')], by = 'geno') |>
            filter(set == fold) |>
            left_join(blue[, c('geno', 'blue')], by = 'geno')

          cvstats[[fold]] = data.frame(
            repl = unique(df2$repl),
            corr = cor(pred$predicted.value, pred$blue, use = 'na.or.complete'),
            mspe = mean((pred$predicted.value - pred$blue)^2, na.rm = TRUE)
          )

        }
        do.call(rbind, cvstats)
      })

      cv.cval[[k]] = cbind(do.call(rbind, lapply(cv, function(x) {
        x$fold = seq_along(x$corr)
        x
      })), c.val = c(c.val))
    }

    cv.cval = do.call(rbind, cv.cval)
    aux = tapply(cv.cval$corr, cv.cval$c.val, mean)
    cvalsel = as.numeric(names(aux[which.max(aux)]))
    
    ablup = ablup |>
      mutate(WABLUP_adj = (1 - c(h2.GA[which(h2.GA$herita == "h2"), 'Estimate'])) /
               (cvalsel + (1 - reli_adj) /
                  reli_adj) * c(h2.GA[which(h2.GA$herita == "h2"), 'Estimate']))
    
    list(
      entries = list(
        BLUE = blue,
        BLUP = blup,
        ABLUP = ablup
      ),
      vc = list(BLUP = vc.GI, ABLUP = vc.GA),
      herita = list(h2.GI, h2.GA),
      lrtest = list(lrtest.GI, lrtest.GA)
    )
  }
)

save(results, file = "Saves/RES_1st.RDA")
