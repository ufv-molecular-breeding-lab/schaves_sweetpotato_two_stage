rm(list=ls())

# Loading the preparations ------------------------------------------------
source("Codes/Preparations.R")
load(file = "Saves/RES_1st.RDA")

# Preparations for each scenario ------------------------------------------
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


# Second-stage models -----------------------------------------------------

## Using BLUEs as entries
blue = transform(blue, geno = as.factor(geno), env = as.factor(env))
Gblue = ASRgenomics::full2sparse(Gmat[rownames(Gmat) %in% levels(blue$geno), colnames(Gmat) %in% levels(blue$geno)])
attr(Gblue, 'INVERSE') = FALSE

mod = asreml(
  fixed = blue ~ env,
  random = ~ vm(geno, Gblue):corh(env) + ide(geno),
  family = asr_gaussian(dispersion = 1),
  weights = "WBLUE",
  data = blue,
  workspace = '1.5gb',
  maxit = 50
)

vc.blue = summary(mod)$varcomp |> mutate(entry = "BLUE", trait = 'rytha') |>
  rownames_to_column('effect')
blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
prblue = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "BLUE", trait = 'rytha')

## Using dBLUPs as entries
blup = blup[-which(is.na(blup$WBLUP)), ]
blup = transform(blup, geno = as.factor(geno), env = as.factor(env))
Gblup = ASRgenomics::full2sparse(Gmat[rownames(Gmat) %in% levels(blup$geno), colnames(Gmat) %in% levels(blup$geno)])
attr(Gblup, 'INVERSE') = FALSE

mod = asreml(
  fixed = dblup ~ env,
  random = ~ vm(geno, Gblup):corh(env) + ide(geno),
  family = asr_gaussian(dispersion = 1),
  weights = "WBLUP",
  data = blup,
  workspace = '1.5gb',
  maxit = 50
)

vc.dblup = summary(mod)$varcomp |> mutate(entry = "dBLUP", trait = 'rytha') |>
  rownames_to_column('effect')
blu = as.data.frame(summary(mod, coef = TRUE)$coef.random)
prdblup = as.data.frame(blu[grep("vm", rownames(blu)), ]) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "dBLUP", trait = 'rytha')

## Using dABLUPs as entries (regular deregression)
ablup = transform(ablup, geno = as.factor(geno), env = as.factor(env))
Gblup = ASRgenomics::full2sparse(Gmat[rownames(Gmat) %in% levels(ablup$geno), colnames(Gmat) %in% levels(ablup$geno)])
attr(Gblup, 'INVERSE') = FALSE

mod = asreml(
  fixed = dAblup ~ env,
  random = ~ vm(geno, Gblup):corh(env),
  family = asr_gaussian(dispersion = 1),
  weights = "WABLUP",
  data = ablup,
  workspace = '1.5gb',
  maxit = 50
)

vc.dablup = summary(mod)$varcomp |> mutate(entry = "dABLUP", trait = 'rytha') |>
  rownames_to_column('effect')
prdablup = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "dABLUP", trait = 'rytha')

## Using dABLUPs as entries (Garrick's et al. deregression)
ablup2 = droplevels(ablup[-which(is.na(ablup$WABLUP_adj)), ])
ablup2 = transform(ablup2, geno = as.factor(geno), env = as.factor(env))
Gblup = ASRgenomics::full2sparse(Gmat[rownames(Gmat) %in% levels(ablup2$geno), colnames(Gmat) %in% levels(ablup2$geno)])
attr(Gblup, 'INVERSE') = FALSE

mod = asreml(
  fixed = dAblup_adj ~ env,
  random = ~ vm(geno, Gblup):corh(env),
  family = asr_gaussian(dispersion = 1),
  weights = "WABLUP_adj",
  data = ablup2,
  workspace = '1.5gb',
  maxit = 50
)
vc.dablup_adj = summary(mod)$varcomp |>
  mutate(entry = "dABLUP_adj", trait = 'rytha') |>
  rownames_to_column('effect')
prdablup_adj = as.data.frame(summary(mod, coef = TRUE)$coef.random) %>%
  rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at("geno", str_replace, ".*_", "") %>%
  mutate(entry = "dABLUP_adj", trait = 'rytha')


# Wrap-up -----------------------------------------------------------------
vc.2step = list(
  BLUE = vc.blue,
  dBLUP = vc.dblup,
  dABLUP = vc.dablup,
  dABLUP_adj = vc.dablup_adj
)
blup.2step = list(
  BLUE = prblue,
  dBLUP = prdblup,
  dABLUP = prdablup,
  dABLUP_adj = prdablup_adj
)

save(vc.2step, blup.2step, file = "Saves/stg2mods.RDA")