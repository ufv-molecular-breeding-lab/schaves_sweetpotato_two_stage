rm(list=ls())

# Loading the preparations ------------------------------------------------
source("Codes/Preparations.R")

# Model ------------------------------------------------------------------
x = pheno[order(pheno$env, pheno$row_number, pheno$col_number), ]
x = x[which(grepl("OT", x$env)),]
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
  check.cod = as.factor(ifelse(status == "Check", geno, NA)),
  env = as.factor(env)
)

ginv = ASRgenomics::G.inverse(G = Gmat, sparseform = FALSE)$Ginv
Hsparse <- H.inverse(A = Amat[which(rownames(Amat) %in% unique(x$geno)),
                              which(colnames(Amat) %in% unique(x$geno))],
                     Ginv = ginv[which(rownames(ginv) %in% unique(x$geno)),
                                 which(colnames(ginv) %in% unique(x$geno))], 
                     sparseform = TRUE, omega = .1, tau = 1- .1)

x$genod = x$geno
mod = asreml(
  fixed = rytha  ~ pool + env + 
    check.cod + check.cod:env,
  random = ~ vm(geno.cod, Hsparse):corh(env) +
    genod:corh(env) +
    rn:at(env) + cn:at(env) + colgroup:at(env) +
    rowgroup:at(env) + rc.group:at(env),
  residual = ~ dsum( ~ ar1v(rn):ar1(cn) | env),
  na.action = na.method(x = "include", y = "include"),
  data = x,
  maxit = 50,
  workspace = '3gb'
)

vc = summary(mod)$varcomp |> mutate(entry = "SS", trait = 'rytha') |>
  rownames_to_column('effect')
blup = as.data.frame(summary(mod, coef = TRUE)$coef.random)[which(grepl("geno.cod", rownames(summary(mod, coef = TRUE)$coef.random))), ] %>% rownames_to_column("geno") %>% select(-z.ratio, -std.error) %>%
  separate("geno", into = c("geno", "env"), sep = ":env_") %>%
  mutate_at('geno', str_replace, ".*\\_", "") %>%
  mutate(entry = "SS", trait = 'rytha')

save(vc, blup, file = "Saves/ssmods.RDA")