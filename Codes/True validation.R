rm(list=ls())

# Loading the preparations ------------------------------------------------
source("Codes/Preparations.R")
load(file = "Saves/RES_1st.RDA")

# Models ------------------------------------------------------------------
x = pheno[which(grepl("PT", pheno$env)),]
# x$rytha[which(x$rytha > 100 | x$rytha < 1)] = NA

x = transform(
  x,
  env = as.factor(env),
  geno = as.factor(geno),
  rn = as.factor(row_number),
  cn = as.factor(col_number),
  pool = as.factor(ifelse(pool == "C", NA, pool)),
  check = as.factor(ifelse(status == "Check", "Check", "Candidate")),
  geno.cod = as.factor(ifelse(status != "Check", geno, NA)),
  check.cod = as.factor(ifelse(status == "Check", geno, NA))
)

Gsparse = full2sparse(Gmat[rownames(Gmat) %in% unique(c(as.character(x$geno.cod))),
                           colnames(Gmat) %in% unique(c(as.character(x$geno.cod)))])
attr(Gsparse, "INVERSE") = FALSE

# unique(x$geno.cod)[!unique(x$geno.cod) %in% attr(Gsparse, 'rowNames')]

modpt = asreml(
  fixed = rytha ~ pool + env + check.cod + check.cod:env,
  random = ~ vm(geno.cod, Gsparse):corh(env) + ide(geno.cod),
  residual = ~ dsum( ~ ar1(rn):ar1(cn) | env),
  na.action = na.method(x = "include", y = "include"),
  data = x,
  workspace = "5gb",
  maxit = 50
)

ptmods = list(
  vc = summary(modpt)$varcomp |> mutate(trait = 'rytha') |>
    rownames_to_column('effect'),
  blup = as.data.frame(summary(modpt, coef = TRUE)$coef.random)[which(grepl("vm", rownames(
    summary(modpt, coef = TRUE)$coef.random)
  )),] %>% rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
    separate("geno", into = c("geno", "env"), sep = ":env_") %>%
    mutate_at('geno', str_extract, paste(paste(attr(Gsparse, 'rowNames'), collapse = "|"))) |>
    mutate(trait = 'rytha')
)

# Real validation ---------------------------------------------------------
ablup = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind,  x$entries[which(names(x$entries) == 'ABLUP')])
  rownames(aux) = NULL
  aux
})) |> rownames_to_column('env') |> mutate_at("env", str_replace, "\\..*",'') |> 
  mutate(entry = "ABLUP") |> rename(WABLUP = WBLUP) |> 
  filter(geno %in% rownames(Gmat)) |> droplevels()
ablup[which(ablup$reli < 0.15),'WABLUP'] = NA
ablup[which(ablup$reli_adj < 0.15),c("dAblup_adj",'WABLUP_adj')] = NA

blup = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind,   x$entries[which(names(x$entries) == 'BLUP')])
  rownames(aux) = NULL
  aux
})) |> rownames_to_column('env') |> mutate_at("env", str_replace, "\\..*",'') |> 
  filter(geno %in% rownames(Gmat)) |> droplevels()

blue = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind,   x$entries[which(names(x$entries) == 'BLUE')])
  rownames(aux) = NULL
  droplevels(aux)
})) |> rownames_to_column('env') |> mutate_at("env", str_replace, "\\..*",'') |> 
  filter(geno %in% rownames(Gmat)) |> droplevels()

## Using BLUEs as entries ---------------------------
blue$yNA = blue$blue
blue[which(blue$geno %in% levels(x$geno)), 'yNA'] = NA
blue = transform(blue, geno = as.factor(geno), env = as.factor(env))
Gblue = full2sparse(Gmat[rownames(Gmat) %in% levels(blue$geno),
                         colnames(Gmat) %in% levels(blue$geno)])
attr(Gblue, 'INVERSE') = FALSE

mod = asreml(
  fixed = yNA ~ env,
  random = ~ vm(geno, Gblue):corh(env) + ide(geno),
  family = asr_gaussian(dispersion = 1),
  weights = "WBLUE",
  data = blue,
  workspace = '1.5gb',
  maxit = 100
)
blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
prblue = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "BLUE", trait = 'rytha') |>
  filter(geno %in% unique(x$geno))

## Using dBLUPs as entries ---------------------------
blup$yNA = blup$dblup
blup[which(blup$geno %in% levels(x$geno)), 'yNA'] = NA
blup = blup[-which(is.na(blup$WBLUP)), ]
blup = transform(blup, geno = as.factor(geno), env = as.factor(env))
Gblup = full2sparse(Gmat[rownames(Gmat) %in% levels(blup$geno), 
                         colnames(Gmat) %in% levels(blup$geno)])
attr(Gblup, 'INVERSE') = FALSE

mod = asreml(
  fixed = yNA ~ env,
  random = ~ vm(geno, Gblup):corh(env) + ide(geno),
  family = asr_gaussian(dispersion = 1),
  weights = "WBLUP",
  data = blup,
  workspace = '1.5gb',
  maxit = 50
)
blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
prdblup = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "BLUP", trait = 'rytha') |>
  filter(geno %in% unique(x$geno))

## Using dABLUPs (regular deregression) as entries ---------------------------
ablup$yNA = ablup$dAblup
ablup[which(ablup$geno %in% levels(x$geno)), 'yNA'] = NA
ablup = transform(ablup, geno = as.factor(geno), env = as.factor(env))
Gblup = full2sparse(Gmat[rownames(Gmat) %in% levels(ablup$geno), 
                         colnames(Gmat) %in% levels(ablup$geno)])
attr(Gblup, 'INVERSE') = FALSE

mod = asreml(
  fixed = yNA ~ env,
  random = ~ vm(geno, Gblup):corh(env),
  family = asr_gaussian(dispersion = 1),
  weights = "WABLUP",
  # na.action = na.method(x = "include", y = "include"),
  data = ablup,
  workspace = '1.5gb',
  maxit = 50
)

prdablup = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "dABLUP", trait = 'rytha')|>
  filter(geno %in% unique(x$geno))

## Using dABLUPs as entries (Garrick's et al. deregression) as entries ---------------------------
ablup2 = droplevels(ablup[-which(is.na(ablup$WABLUP_adj)), ])
ablup2$yNA = ablup2$dAblup_adj
ablup2[which(ablup2$geno %in% levels(x$geno)), 'yNA'] = NA
ablup2 = transform(ablup2, geno = as.factor(geno), env = as.factor(env))
Gblup = full2sparse(Gmat[rownames(Gmat) %in% levels(ablup2$geno), colnames(Gmat) %in% levels(ablup2$geno)])
attr(Gblup, 'INVERSE') = FALSE

mod = asreml(
  fixed = yNA ~ env,
  random = ~ vm(geno, Gblup):corh(env),
  family = asr_gaussian(dispersion = 1),
  weights = "WABLUP_adj",
  # na.action = na.method(x = "include", y = "include"),
  data = ablup2,
  workspace = '1.5gb',
  maxit = 50
)

prdablup_adj = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "dABLUP_adj", trait = 'rytha') |>
  filter(geno %in% unique(x$geno))

blup.rval = list(BLUE = prblue, dBLUP = prdblup, dABLUP = prdablup,
                 dABLUP_adj = prdablup_adj)

save(blup.rval, file = 'Saves/rval_pred.RDA')


