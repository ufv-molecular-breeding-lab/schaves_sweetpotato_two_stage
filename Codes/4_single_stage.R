source("Codes/1_preparations.R")
source("Codes/Functions/measure_variances.R")
source("Codes/Functions/expvar.R")
source("Codes/Functions/fa_summa_rr.R")
library(ComplexHeatmap)

# Further preparations ------------
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

# Single-stage model ------------------------------------------------------
famod_sing = list()
for (i in 1:4) {
  a = Sys.time()
  famod_sing[[paste0("FA", i)]] = asreml(
    fixed = rytha  ~ pool + env*check.cod,
    random = ~ rr(env, i):vm(genokeep, ginv) + diag(env):vm(genokeep, ginv) +
      corh(env):ide(genokeep) +
      at(env):colgroup + at(env):rowgroup + at(env):colgroup:rowgroup,
    residual = ~ dsum(~ar1v(rn):ar1(cn)|env) ,
    sparse = ~ env:genodrop,
    na.action = na.method(x = "include", y = "include"),
    data = dat,
    maxit = 30,
    workspace = '5gb'
  )
  if(any(na.exclude(famod_sing[[paste0("FA", i)]]$vparameters.pc > 1))) famod_sing[[paste0("FA", i)]] = up.mod(famod_sing[[paste0("FA", i)]])
  b = Sys.time()
  famod_sing[[paste0("FA", i)]]$ellapsed = b-a
  if(!famod_sing[[paste0("FA", i)]]$converge) break
}
save(famod_sing, file = "Saves/famod_sing.RDA")

# Model selection ---------------------------------------------------------
load(file = "Saves/famod_sing.RDA")
modsel = lapply(famod_sing, expvar, env = 'env', geno = 'genokeep')

temp = data.frame(
  model = seq_along(famod_sing),
  ellapsed = do.call(c, lapply(famod_sing, function(x) x$ellapsed))
)

fasingsel = famod_sing[[2]]
rm(famod_sing)

# Working on the selected model -------------------------------------------
name.env = levels(fasingsel$mf$env)
num.env = nlevels(fasingsel$mf$env)
name.geno = levels(fasingsel$mf$genokeep)
num.geno = nlevels(fasingsel$mf$genokeep)
summa = fa.summa.rr(model = fasingsel, env = 'env', geno = 'genokeep')
cor_hsc = summa$nav[grep('cor', rownames(summa$nav)),1]
s2na = summa$nav[!grepl('cor', rownames(summa$nav)),1]
G_na = sqrt(diag(s2na)) %*% 
  (diag(length(s2na)) + 
     cor_hsc*(matrix(1, nrow = length(s2na), ncol = length(s2na)) - diag(length(s2na)))) %*% 
  sqrt(diag(s2na)) 

measure_variances(summa$Gvcov)

# Heatmap with genetic correlations ----------------------------------------------------
colfun = circlize::colorRamp2(breaks = c(-.2, 0, 1), colors = c("firebrick", "white", "forestgreen"), 
                    space = "RGB")
draw(
  Heatmap(
    summa$Gcor, 
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    show_heatmap_legend = FALSE,
    row_split = case_when(
      grepl("NAM", rownames(summa$Gcor)) ~ "Namulonge",
      grepl("SER", rownames(summa$Gcor)) ~ "Serere",
      grepl("RWE", rownames(summa$Gcor)) ~ "Rwebitaba"
    ),
    column_split = case_when(
      grepl("NAM", colnames(summa$Gcor)) ~ "Namulonge",
      grepl("SER", colnames(summa$Gcor)) ~ "Serere",
      grepl("RWE", colnames(summa$Gcor)) ~ "Rwebitaba"
    ),
    col = colfun, 
    column_labels = str_extract(colnames(summa$Gcor), "22|21"), 
    row_labels = str_extract(rownames(summa$Gcor), "22|21"), 
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.3f", summa$Gcor[i, j]), x, y, gp = gpar(fontsize = 10))
    }
  ),
  heatmap_legend_side = "top"
)
