rm(list=ls())

# Preps -----------------------------------------------------------------
source("Codes/Preparations.R")
load(file = "Saves/RES_1st.RDA")
library(ggrepel)
library(ggpubr)
library(GGally)

# Phenotype: Overview ----------------------------------------------------------------
pheno |> 
  mutate(env2 = env) |> 
  separate('env2', into = c('loc', 'year'), sep = 6) |> 
  mutate(loc = gsub("OT[0-9]|PT[0-9]",'', loc), 
         year = paste0("20", year),
         # year = ifelse(grepl("PT3", env), 2024, year),
         env = factor(env, levels = c("OT1RWE21", "OT2NAM21", "OT2SER21",
                                      "OT1RWE22", "OT1NAM22", "OT1SER22",
                                      "PT1RWE23","PT1NAM23", "PT1SER23",
                                      "PT3RWE23", "PT3NAM23", "PT3SER23"))) |> 
  ggplot(aes(y = rytha, x = env)) + 
  geom_boxplot(aes(fill = loc)) +
  facet_grid(~year, scale = "free") + 
  labs(x = "Environment", y = "Storage root yield (ton/ha)") + 
  theme_bw() + 
  theme(legend.title = element_blank(), legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) + 
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada'))

# Phenotype: Grid --------------------------------------------------------------------
pheno$status2 = pheno$status
pheno[which(is.na(pheno$rytha)), "status2"] = NA
ggplot() +
  geom_tile(
    data = pheno[which(pheno$status == "Check" &
                         grepl("OT", pheno$env)), ],
    aes(x = col_number, y = row_number, fill = ifelse(is.na(rytha), NA, geno)),
    # color = "black",
    size = 0.5,
    width = .9,
    height = .9
  ) +
  geom_tile(
    data = pheno[which(pheno$status != "Check" &
                         grepl("OT", pheno$env)), ],
    aes(x = col_number, y = row_number, fill = status2),
    # color = "black",
    size = 0.5,
    width = .9,
    height = .9
  )  +
  facet_wrap(. ~ factor(
    env,
    levels = c(
      "OT2NAM21",
      "OT2SER21",
      "OT1RWE21",
      "OT1NAM22",
      "OT1SER22",
      "OT1RWE22"
    )
  )) +
  scale_fill_manual(
    breaks = c("Parent", "Regular", unique(pheno[which(pheno$status == "Check" &
                                                         grepl("OT", pheno$env)), 'geno'])),
    values = c(
      "Parent" = 'red3',
      "Regular" = 'grey',
      "Tanzania" = '#377eb8',
      "New_Kawogo" = '#4daf4a',
      "Cemsa_74-228" = '#984ea3',
      "NAROSPOT1" = '#ff7f00',
      "Beauregard" = '#ffff33',
      "NASPOT_8" = '#a65628',
      "Silk_Omuyaka" = '#f781bf',
      "RAK808" = "black"
    ),
    na.value = "white"
  ) +
  labs(y = "Row", x = "Column", fill = "Genotype") +
  theme_bw() + theme(legend.position = 'right') +
  layer(
    geom = "segment",
    stat = StatTileBorder,
    data = pheno[which(grepl("OT", pheno$env)), ],
    mapping = aes(
      group = 1,
      grp = rowgroup,
      y = row_number,
      x = col_number
    ),
    position = 'identity'
  ) +
  layer(
    geom = "segment",
    stat = StatTileBorder,
    data = pheno[which(grepl("OT", pheno$env)), ],
    mapping = aes(
      group = 1,
      grp = colgroup,
      y = row_number,
      x = col_number
    ),
    position = 'identity'
  )

temp = pheno[which(grepl("PT", pheno$env)), ]
temp = left_join(
  temp,
  do.call(rbind, lapply(split(temp, temp$env), function(x)
    as.data.frame(table(x$geno)))) |>
    rownames_to_column('env') |> mutate_at('env', str_replace, '\\..*', ''),
  by = c("env", 'geno' = "Var1")
)

ggplot() +
  geom_tile(
    data = temp[which(temp$status == "Check"), ],
    aes(
      x = col_number,
      y = row_number,
      fill = geno,
      alpha = factor(Freq)
    ),
    color = "black",
    size = 0.5,
    width = .9,
    height = .9
  ) +
  geom_tile(
    data = temp[which(temp$status != "Check"), ],
    aes(
      x = col_number,
      y = row_number,
      fill = status2,
      alpha = factor(Freq)
    ),
    color = "black",
    size = 0.5,
    width = .9,
    height = .9
  ) +
  scale_fill_manual(
    breaks = c("Parent", "Regular", unique(temp[which(temp$status == "Check"), 'geno'])),
    values = c(
      "Parent" = 'red3',
      "Regular" = 'black',
      "Tanzania" = '#377eb8',
      "New_Kawogo" = '#4daf4a',
      "Cemsa_74-228" = '#984ea3',
      "NAROSPOT1" = '#ff7f00',
      "Beauregard" = '#ffff33',
      "NASPOT_8" = '#a65628',
      "Silk_Omuyaka" = '#f781bf'
    ),
    na.value = "white"
  ) +
  scale_alpha_manual(values = c("1" = 0.3, "2" = 0.6, "3" = 1),
                     breaks = as.character(unique(temp$Freq))) +
  theme_bw() + theme(legend.position = 'right') +
  scale_x_discrete(limits = as.character(seq(1:11))) +
  scale_y_discrete(limits = as.character(seq(1:15))) +
  labs(x = "Columns",
       y = "Rows",
       fill = "Genotype",
       alpha = "No. replicates") +
  facet_wrap(. ~ env)
rm(temp)

# G matrix: PCA -----------------------------------------------------------------
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

PC = FactoMineR::PCA(X = Gmat, graph = FALSE, scale.unit = FALSE)
pc.plot = data.frame(
  geno = names(PC$ind$coord[,1]),
  pc1 = PC$ind$coord[,1],
  pc2 = PC$ind$coord[,2], row.names = NULL
) |> left_join(unique(pheno[,c("geno",'pool','status')]), by = "geno")


pca.A = FactoMineR::PCA(
  Gmat[
    rownames(Gmat) %in% c(unique(pheno[which(pheno$pool == "A"),'geno']),
                          unique(pheno[which(pheno$status == "Parent" & pheno$pool == "A"), 'geno'])),
    colnames(Gmat) %in% c(unique(pheno[which(pheno$pool == "A"),'geno']),
                          unique(pheno[which(pheno$status == "Parent" & pheno$pool == "A"), 'geno']))
  ], scale.unit = TRUE, graph = FALSE)


pca.B = FactoMineR::PCA(
  Gmat[
    rownames(Gmat) %in% c(unique(pheno[which(pheno$pool == "B"),'geno']),
                          unique(pheno[which(pheno$status == "Parent" & pheno$pool == "B"), 'geno'])),
    colnames(Gmat) %in% c(unique(pheno[which(pheno$pool == "B"),'geno']),
                          unique(pheno[which(pheno$status == "Parent" & pheno$pool == "B"), 'geno']))
  ], scale.unit = TRUE, graph = FALSE)

pc.df.A = data.frame(
  geno = names(pca.A$ind$coord[,1]),
  pc1 = pca.A$ind$coord[,1], 
  pc2 = pca.A$ind$coord[,2], row.names = NULL
) |> left_join(pc.plot[,c('geno', 'status')])

pc.df.B = data.frame(
  geno = names(pca.B$ind$coord[,1]),
  pc1 = pca.B$ind$coord[,1], 
  pc2 = pca.B$ind$coord[,2], row.names = NULL
) |> left_join(pc.plot[,c('geno', 'status')])

ggarrange(
  pc.plot |> filter(status != "Check") |>
    ggplot(aes(x = pc1, y = pc2)) +
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_point(
      aes(color = pool, shape = status),
      size = 2,
      alpha = .8
    ) +
    geom_label_repel(
      data = subset(pc.plot, subset = status == "Parent"),
      aes(label = geno, fill = pool),
      size = 2,
      box.padding = .5,
      alpha = .7,
      fontface = 'bold',
      show.legend = FALSE
    ) +
    labs(
      x = paste0('PC1 (', round(PC$eig[1, 2], 2), '%)'),
      title = "Complete dataset",
      y = paste0('PC2 (', round(PC$eig[2, 2], 2), '%)')
    ) +
    scale_color_manual("Pool", values = c('#1b9e77', '#d95f02')) +
    scale_fill_manual("Pool", values = c('#1b9e77', '#d95f02'))  +
    scale_shape_manual("", values = c(15, 17)) +
    theme_bw() + theme(legend.position = 'bottom'),
  ggarrange(
    ggplot(data = pc.df.A, aes(x = pc1, y = pc2))  +
      geom_point(colour = '#1b9e77', aes(shape = status))  +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(
        x = paste0('PC1 (', round(pca.A$eig[1, 2], 2), '%)'),
        y = paste0('PC2 (', round(pca.A$eig[2, 2], 2), '%)'),
        shape = '',
        alpha = '',
        title = "Pool A"
      ) +
      scale_shape_manual("", values = c(15, 17)) +
      ylim(-30, 30) + xlim(-30, 30) ,
    ggplot(data = pc.df.B, aes(x = pc1, y = pc2))  +
      geom_point(colour = '#d95f02', aes(shape = status))  +
      theme_bw() +
      theme(legend.position = 'none') +
      labs(
        x = paste0('PC1 (', round(pca.B$eig[1, 2], 2), '%)'),
        y = paste0('PC2 (', round(pca.B$eig[2, 2], 2), '%)'),
        shape = '',
        alpha = '',
        title = "Pool B"
      ) + ylim(-40, 25) + xlim(-40,25) +
      scale_shape_manual("", values = c(15, 17))), nrow = 2, common.legend = TRUE
)

# First stage: LRT --------------------------------------------------------
lrtest = do.call(rbind, lapply(results, function(x){
  x$lrtest[[1]] = x$lrtest[[1]] |> mutate(effect = "G")
  colnames(x$lrtest[[1]]) = gsub("\\-|\\(|\\)",".",colnames(x$lrtest[[1]]))
  aux = do.call(rbind, x$lrtest)
  rownames(aux) = NULL
  aux
})) |> rownames_to_column("env") |> mutate_at("env", str_replace, '\\..*','') |> 
  mutate(effect = case_when(effect == "NA" ~ 'N-A', .default = effect))

lrtest %>%
  ggplot(aes(x = env, y = -log(Pr.Chisq.), fill = effect)) + 
  geom_col(position = position_dodge(), colour = "black") + 
  facet_wrap(.~entry) +
  geom_hline(aes(yintercept = -log(0.05), linetype = "p-value = 0.05")) + 
  labs(x = "Environment", y = "-log(P-value) (chi-square test)", linetype = "",
       fill = "Effect") +
  scale_linetype_manual(values = "dashed") + 
  theme_bw() + theme(legend.position = 'top', axis.text.x = element_text(angle = 90)) 

# First stage: Heritability and variance components -----------------------------------------
herita = do.call(rbind, lapply(results, function(x){
  aux = do.call(rbind, x$herita)
  rownames(aux) = NULL
  aux
})) |> rownames_to_column("env") |> mutate_at("env", str_replace, '\\..*','')

vc = do.call(rbind, lapply(results, function(x) do.call(rbind, x$vc))) |> 
  rownames_to_column('env') |> mutate_at('env', str_replace, '\\..*', '') |> 
  filter(!grepl("cor|!R", effect)) |> 
  mutate(effect = factor(case_when(
    effect == 'colgroup' ~ "Col.group",
    effect == "rowgroup" ~ "Row.group",
    effect == 'rn' ~ "Row", 
    effect == 'cn' ~ "Col",
    effect == "genoadd" ~ "Additive",
    effect == "genod" ~ "Non-additive",
    effect == 'geno' ~ "Genotypic",
    effect == 'rc.group' ~ "Block",
    grepl("var", effect) ~ "Error"
  ), levels = c("Genotypic", "Additive", "Non-additive", "Row", "Col", "Row.group",
                "Col.group", "Block", "Error")))

ggarrange(
  ggplot(data = vc |> filter(effect %in% c("Genotypic", "Additive", "Non-additive")),
         aes(x = model, y = component)) + 
    geom_col(aes(fill = effect)) + 
    theme_bw() + 
    facet_wrap(env ~ ., scales = "free") + 
    theme(legend.position = 'top') + 
    labs(x = "Model", y = "Variance component estimate", fill = "Effect") + 
    # geom_errorbar(aes(ymin = component - std.error, ymax = component + std.error),
    #               width = .2) + 
    # geom_text(data = subset(vc, subset = effect %in% c("Genotypic", "Additive", "Non-additive")),
    #           label = '*', aes(y = component + std.error), vjust = .1) +
    scale_fill_manual(values = c('steelblue', 'firebrick', 'forestgreen')),
  herita |> mutate(herita = factor(case_when(herita == "h2" ~ "Narrow-sense",
                                             herita == "H2" ~ "Broad-sense"),
                                   levels = c("Narrow-sense", "Broad-sense"))) |> 
    ggplot(aes(y = Estimate, x = env, fill = herita)) + 
    facet_grid(. ~ entry) + 
    theme_bw() + 
    theme(legend.position = 'top', 
          axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
    geom_col(position = position_dodge(), colour = "black")  +
    scale_fill_manual(values = c("#8da0cb", "#fc8d62")) + 
    labs(x = "Environments", fill = "Heritability") + 
    geom_errorbar(aes(ymax = Estimate + SE, ymin = Estimate - SE),
                  position = position_dodge(.9), width = .3),
  labels = "AUTO", nrow = 1
)

# First stage: reliabilities ----------------------------------------------
relis = rbind(
  ablup[,c("env", 'geno', 'reli', 'reli_adj', 'entry')] |> 
    pivot_longer(reli:reli_adj, names_to = "reli") |> 
    mutate(reli = case_when(reli == "reli_adj" ~ "Adjusted", .default = "Standard")),
  blup[,c('env', 'geno', 'reli', 'entry')] |> mutate(value = reli, reli = "Standard")
)

entry = full_join(ablup[,c("env", 'geno', 'dAblup', 'dAblup_adj', 'Ablup')],
                  blup[,c('env', 'geno', 'dblup', 'blup')], by = c("env", 'geno')) |> 
  full_join(blue[,c("env", 'geno','blue')], by = c("env", 'geno'))

ggplot(data = relis, aes(x = entry, y = value, fill = reli)) + 
  geom_violin() +
  theme_bw() + 
  geom_hline(aes(yintercept = 0.15), linetype = 'dashed') +
  theme(legend.position = 'top') + 
  labs(x = "Entry", y = "Value", fill = "Reliability") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3"))


# First stage: Means -----------------------------------------------------
entry = full_join(ablup[, c("env", 'geno', 'dAblup', 'dAblup_adj', 'Ablup')],
                  blup[, c('env', 'geno', 'dblup', 'blup')], by = c("env", 'geno')) |>
  full_join(blue[, c("env", 'geno', 'blue')], by = c("env", 'geno'))


lowerFn <- function(data, mapping, emap = NULL, ...) {
  mapping <- ggplot2:::new_aes(c(mapping, emap))
  xx = rlang::quo_get_expr(mapping$x)
  yy = rlang::quo_get_expr(mapping$y)
  l1 = quantile(data[,paste(xx)], probs = .90, na.rm = TRUE)
  l2 = quantile(data[,paste(yy)], probs = .90, na.rm = TRUE)
  
  dat1 = data[which(!data[,paste(xx)] >= l1 | !data[,paste(yy)] >= l2),]
  
  dat2 = data[which(data[,paste(xx)] >= l1 & data[,paste(yy)] >= l2),]
  dat3 = data.frame(y = min(dat2[,paste(yy)]), x = min(dat2[,paste(xx)]))
  
  p = ggplot(mapping = mapping) +
    geom_point(data = data,
               colour = 'black', size = 1.5, show.legend = F, alpha = .4) +
    geom_rug(data = data,
             colour = 'black')  + theme_bw()
}
diagFn <- function(data, mapping, emap = NULL, ...) {
  mapping <- ggplot2:::new_aes(c(mapping, emap))
  xx = rlang::quo_get_expr(mapping$x)
  l1 = quantile(data[,paste(xx)], probs = .90, na.rm = TRUE)
  
  dat1 = data[which(data[,paste(xx)] >= l1),]
  dat2 = data[which(!data[,paste(xx)] >= l1),]
  
  p = ggplot(mapping = mapping) + 
    geom_density(data = data, colour = "black", fill = "black", alpha = .3) +
    theme_bw()
}
upperFn <- function(data, mapping, size = 3, ...) {
  
  # get the x and y data to use the other code
  xx = rlang::quo_get_expr(mapping$x)
  yy = rlang::quo_get_expr(mapping$y)
  
  ct <- cor(data[,paste(xx)], data[,paste(yy)], method = "spearman", 
            use = "na.or.complete")
  rt <- format(ct, digits=3)
  
  l1 = quantile(data[,paste(xx)], probs = .90, na.rm = TRUE)
  l2 = quantile(data[,paste(yy)], probs = .90, na.rm = TRUE)
  
  dat1 = data[which(!data[,paste(xx)] >= l1 | !data[,paste(yy)] >= l2),]
  
  dat2 = data[which(data[,paste(xx)] >= l1 & data[,paste(yy)] >= l2),]
  dat3 = data.frame(y = min(dat2[,paste(yy)]), x = min(dat2[,paste(xx)]))
  
  coinc = paste0("Coinc.: ",round(nrow(dat2)/(nrow(data) * 0.1) * 100, 2), '%')
  
  # plot the cor value
  ggally_text(
    label = paste("Rank corr.:", rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.55, 
    color = "black",
    size = size,
    ...
  ) +
    theme_classic() + 
    theme(
      panel.background = element_rect(
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}

ggpairs(
  data = entry[, -c(1, 2)] |>
    rename(
      dBLUP = dblup,
      BLUE = blue,
      dABLUP = dAblup,
      `dABLUP*` = dAblup_adj,
      BLUP = blup,
      ABLUP = Ablup
    ) |>
    dplyr::select(BLUE, BLUP, ABLUP, dBLUP, dABLUP, `dABLUP*`),
  lower = list(continuous = wrap(lowerFn)),
  diag = list(continuous = wrap(diagFn)),
  upper = list(continuous = wrap(upperFn))
) + theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))


# Second stage and single-stage: Entries ----------------------------------
load(file = "Saves/ssmods.RDA")
load(file = "Saves/stg2mods.RDA")

lowerFn <- function(data, mapping, emap=NULL, int = .1, ...) {
  mapping <- ggplot2:::new_aes(c(mapping, emap))
  xx = rlang::quo_get_expr(mapping$x)
  yy = rlang::quo_get_expr(mapping$y)
  top.x = data[order(data[,paste(xx)], decreasing = TRUE), 1][1:(ceiling(nrow(data)*int))]
  top.y =  data[order(data[,paste(yy)], decreasing = TRUE), 1][1:(ceiling(nrow(data)*int))]
  top.xy = intersect(top.x, top.y)
  minmu.x = data[order(data[,paste(xx)], decreasing = TRUE), paste(xx)][ceiling(nrow(data)*int)]
  minmu.y = data[order(data[,paste(yy)], decreasing = TRUE), paste(yy)][ceiling(nrow(data)*int)]
  
  dat1 = data[which(!data[,1] %in% top.xy),]
  
  dat2 = data[which(data[,1] %in% top.xy),]
  dat3 = data.frame(y = minmu.y, x = minmu.x)
  
  p = ggplot(mapping = mapping) + 
    geom_point(data = dat1,
               colour = 'black', size = 1.7, show.legend = F, alpha = .3) +
    geom_point(data = dat2,
               colour = '#0075cc', size = 1.7, show.legend = F, alpha = .4) +
    geom_rug(data = dat1,
             colour = 'black') +
    geom_rug(data = dat2,
             colour = '#0075cc')
}
diagFn <- function(data, mapping, emap=NULL, ...) {
  mapping <- ggplot2:::new_aes(c(mapping, emap))
  xx = rlang::quo_get_expr(mapping$x)
  l1 = quantile(data[,paste(xx)], probs = .90, na.rm = TRUE)
  
  dat1 = data[which(data[,paste(xx)] >= l1),]
  dat2 = data[which(!data[,paste(xx)] >= l1),]
  
  p = ggplot(mapping = mapping) + 
    geom_density(data = dat2, colour = "black", fill = "black", alpha = .3) +
    geom_density(data = dat1, colour = "#0075cc", fill = "#0075cc", alpha = .3) + 
    theme_bw()
}
upperFn <- function(data, mapping, size = 3, int = .1, ...) {
  
  # get the x and y data to use the other code
  xx = rlang::quo_get_expr(mapping$x)
  yy = rlang::quo_get_expr(mapping$y)
  
  ct <- cor(data[,paste(xx)], data[,paste(yy)], method = "spearman", 
            use = "complete.obs")
  rt <- format(ct, digits=3)
  
  top.x = data[order(data[,paste(xx)], decreasing = TRUE), 1][1:(ceiling(nrow(data)*int))]
  top.y =  data[order(data[,paste(yy)], decreasing = TRUE), 1][1:(ceiling(nrow(data)*int))]
  top.xy = intersect(top.x, top.y)
  minmu.x = data[order(data[,paste(xx)], decreasing = TRUE), paste(xx)][ceiling(nrow(data)*int)]
  minmu.y = data[order(data[,paste(yy)], decreasing = TRUE), paste(yy)][ceiling(nrow(data)*int)]
  
  dat1 = data[which(!data[,1] %in% top.xy),]
  
  dat2 = data[which(data[,1] %in% top.xy),]
  dat3 = data.frame(y = minmu.y, x = minmu.x)
  
  coinc = paste0("Coinc.: ", round(length(top.xy)/ceiling(nrow(data)*int) * 100, 2), '%')
  
  # plot the cor value
  ggally_text(
    label = paste("Rank corr.:", rt), 
    mapping = aes(),
    xP = 0.5, yP = 0.55, 
    color = "black",
    size = size,
    ...
  ) + 
    # Coincidence
    geom_text(
      aes(
        x = 0.503,
        y = 0.4
      ),
      label = coinc, 
      color = '#0075cc', fontface = 'bold',
      size = size,
      ...
    ) + 
    # remove all the background stuff and wrap it with a dashed line
    theme_classic() + 
    theme(
      panel.background = element_rect(
        linetype = "longdash"
      ), 
      axis.line = element_blank(), 
      axis.ticks = element_blank(), 
      axis.text.y = element_blank(), 
      axis.text.x = element_blank()
    )
}

met.blup = rbind(
  do.call(rbind, blup.2step),
  blup 
)

GGally::ggpairs(data = met.blup |> reframe(gebv = mean(solution), 
                                           .by = c('geno','entry')) |> 
                  mutate(entry = case_when(
                    entry == "dABLUP_adj" ~ "2S-dABLUP*",
                    entry == "BLUE" ~ "2S-BLUE",
                    entry == 'dBLUP' ~ "2S-dBLUP",
                    entry == "dABLUP" ~ '2S-dABLUP',
                    entry == "SS" ~ "1S"
                  ))|> 
                  pivot_wider(values_from = gebv, names_from = entry), 
                columns = 2:6,
                lower = list(continuous = wrap(lowerFn)),
                diag = list(continuous = wrap(diagFn)),
                upper = list(continuous = wrap(upperFn)))

# Cross-validations across environments ------------------------------
load('Saves/cv2_intrapool.RDA')
load('Saves/cv1_intrapool.RDA')

cvres.intra = rbind(
  do.call(rbind, lapply(cv2$A, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$set), function(x) {
        data.frame(
          corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
          mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
        )
      })) |> rownames_to_column("set")
    })) |> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "A",
      cv = "CV2"
    ),
  do.call(rbind, lapply(cv2$B, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$set), function(x) {
        data.frame(
          corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
          mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
        )
      })) |> rownames_to_column("set")
    })) |> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "B",
      cv = "CV2"
    ),
  do.call(rbind, lapply(cv1$A, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$set), function(x) {
        data.frame(
          corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
          mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
        )
      })) |> rownames_to_column("set")
    })) |> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "A",
      cv = "CV1"
    ),
  do.call(rbind, lapply(cv1$B, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$set), function(x) {
        data.frame(
          corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
          mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
        )
      })) |> rownames_to_column("set")
    })) |> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "B",
      cv = "CV1"
    )
)

load('Saves/cv2_twopool.RDA')
load('Saves/cv1_twopool.RDA')

cvres.inter = rbind(
  do.call(rbind, lapply(cv2, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$set), function(x) {
        data.frame(
          corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
          mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
        )
      })) |> rownames_to_column("set")
    })) |> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      cv = "CV2"
    ),
  do.call(rbind, lapply(cv1, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$set), function(x) {
        data.frame(
          corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
          mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
        )
      })) |> rownames_to_column("set")
    })) |> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      cv = "CV1"
    )
)

col_label = c("Pool A", "Pool B")
names(col_label) = c("A", "B")

entrylab = c("BLUE" = "2S-BLUE", 'dABLUP_adj' = '2S-dABLUP*', 
             'dABLUP' = '2S-dABLUP', 'dBLUP' = '2S-dBLUP')
entrycol = c('BLUE' = '#7fc97f', 'dABLUP_adj' = '#beaed4',
             'dABLUP' = '#fdc086', 'dBLUP' = '#386cb0')

ggarrange(
  ggarrange(
    cvres.inter |>  #reframe(corr = mean(corr), mspe = mean(mspe), .by = c(entry, cv, repl)) |> 
      mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj"))) |> #reframe(corr = mean(corr), .by = c("entry",'cv','repl')) |> 
      ggplot(aes(x = entry, y = corr, fill = entry)) + 
      facet_wrap(cv~.) + 
      geom_boxplot() +
      theme_bw() + theme(legend.position = "top", axis.text.x = element_blank(), 
                         axis.title.x = element_blank())  +
      labs(x = "Pool", y = "Predictive ability", fill = "")  +
      scale_fill_manual(values = entrycol, labels = entrylab) +
      guides(fill=guide_legend(nrow=1)) + 
      ylim(0.05, 0.39),
    
    cvres.inter |>  #reframe(corr = mean(corr), mspe = mean(mspe), .by = c(entry, cv, repl)) |> 
      mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj"))) |> #reframe(mspe = mean(mspe), .by = c("entry",'cv','repl')) |> 
      ggplot(aes(x = entry, y = mspe, fill = entry)) + 
      facet_wrap(cv~.) + 
      geom_boxplot() +
      theme_bw() + theme(legend.position = "top", axis.text.x = element_blank(), 
                         axis.title.x = element_blank()) +
      labs(x = "Pool", y = "MSPE", fill = "")  +
      scale_fill_manual(values = entrycol, labels = entrylab) +
      guides(fill=guide_legend(nrow=1)) + 
      ylim(290, 450),
    common.legend = TRUE,  nrow = 1
  ),
  ggarrange(
    cvres.intra |> #reframe(corr = mean(corr), mspe = mean(mspe), .by = c(entry, pool, cv, repl)) |> 
      mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj"))) |>  
      ggplot(aes(x = entry, y = corr, fill = entry)) + 
      facet_grid(pool~cv,
                 labeller = labeller(.rows = col_label)) + 
      geom_boxplot() +
      theme_bw() + theme(legend.position = "none", axis.text.x = element_blank(),
                         axis.title.x = element_blank()) +
      labs(y = "Predictive ability", fill = "") +
      scale_fill_manual(values = entrycol, labels = entrylab),
    
    cvres.intra |> #reframe(corr = mean(corr), mspe = mean(mspe), .by = c(entry, pool, cv, repl)) |> 
      mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj"))) |>  
      #reframe(mspe = mean(mspe), .by = c("entry",'pool','cv','repl')) |> 
      ggplot(aes(x = entry, y = mspe, fill = entry)) + 
      facet_grid(pool~cv, scales = "free_y",
                 labeller = labeller(.rows = col_label)) + 
      geom_boxplot() +
      theme_bw() + theme(legend.position = "none", axis.text.x = element_blank(),
                         axis.title.x = element_blank()) +
      labs(y = "MSPE", fill = "") +
      scale_fill_manual(values = entrycol, labels = entrylab),  nrow = 1
  ),
  nrow = 2, common.legend = TRUE, labels = "AUTO"
)

# Cross-validations within environments -----------------------------------
rm(list=ls()[-which(grepl("entry", ls()))])
load('Saves/cv2_intrapool.RDA')
load('Saves/cv1_intrapool.RDA')

cvres.intra = rbind(
  do.call(rbind, lapply(cv2$A, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$env), function(x){
        do.call(rbind, lapply(split(x, x$set), function(x) {
          data.frame(
            corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
            mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
          )
        })) |> rownames_to_column("set")
      })) |> rownames_to_column("env") |> mutate_at('env', str_replace, '\\..','')
    })) #|> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "A",
      cv = "CV2"
    ),
  do.call(rbind, lapply(cv2$B, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$env), function(x){
        do.call(rbind, lapply(split(x, x$set), function(x) {
          data.frame(
            corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
            mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
          )
        })) |> rownames_to_column("set")
      })) |> rownames_to_column("env") |> mutate_at('env', str_replace, '\\..','')
    })) #|> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "B",
      cv = "CV2"
    ),
  do.call(rbind, lapply(cv1$A, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$env), function(x){
        do.call(rbind, lapply(split(x, x$set), function(x) {
          data.frame(
            corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
            mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
          )
        })) |> rownames_to_column("set")
      })) |> rownames_to_column("env") |> mutate_at('env', str_replace, '\\..','')
    })) #|> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "A",
      cv = "CV1"
    ),
  do.call(rbind, lapply(cv1$B, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$env), function(x){
        do.call(rbind, lapply(split(x, x$set), function(x) {
          data.frame(
            corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
            mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
          )
        })) |> rownames_to_column("set")
      })) |> rownames_to_column("env") |> mutate_at('env', str_replace, '\\..','')
    })) #|> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      pool = "B",
      cv = "CV1"
    )
)

ggarrange(
  cvres.intra |>
    mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj"))) |> #reframe(corr = mean(corr), .by = c("entry",'pool','cv','repl')) |> 
    ggplot(aes(x = pool, y = corr, fill = entry)) + 
    facet_grid(env ~ cv, scales = "free_y") + 
    geom_boxplot(outlier.size = .5, linewidth = .4) +
    theme_bw() + theme(legend.position = "top") +
    labs(x = "Pool", y = "Predictive ability", fill = "") +
    scale_fill_manual(values = entrycol, labels = entrylab),
  
  cvres.intra |>
    mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj"))) |> #reframe(corr = mean(corr), .by = c("entry",'pool','cv','repl')) |> 
    ggplot(aes(x = pool, y = mspe, fill = entry)) + 
    facet_grid(env~ cv, scales = "free_y") + 
    geom_boxplot(outlier.size = .5, linewidth = .4) +
    theme_bw() + theme(legend.position = "top") +
    labs(x = "Pool", y = "MSPE", fill = "")  +
    scale_fill_manual(values = entrycol, labels = entrylab),
  common.legend = TRUE, labels = "AUTO", nrow = 1
)

rm(list=ls()[-which(grepl("entry", ls()))])
load('Saves/cv2_twopool.RDA')
load('Saves/cv1_twopool.RDA')

cvres.inter = rbind(
  do.call(rbind, lapply(cv2, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$env), function(x){
        do.call(rbind, lapply(split(x, x$set), function(x) {
          data.frame(
            corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
            mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
          )
        })) |> rownames_to_column("set")
      })) |> rownames_to_column("env") |> mutate_at('env', str_replace, '\\..','')
    })) #|> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      cv = "CV2"
    ),
  do.call(rbind, lapply(cv1, function(x) {
    do.call(rbind, lapply(x, function(x) {
      do.call(rbind, lapply(split(x, x$env), function(x){
        do.call(rbind, lapply(split(x, x$set), function(x) {
          data.frame(
            corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
            mspe = mean((x$ybar - x$gebv) ^ 2, na.rm = TRUE)
          )
        })) |> rownames_to_column("set")
      })) |> rownames_to_column("env") |> mutate_at('env', str_replace, '\\..','')
    })) #|> mutate(repl = rep(1:5, each = 5))
  })) |> rownames_to_column('entry') |>
    mutate(
      entry = str_replace(entry, '\\..*', ''),
      cv = "CV1"
    )
) 

ggarrange(
  cvres.inter |>
    mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj")))|> #reframe(corr = mean(corr), .by = c("entry",'pool','cv','repl')) |> 
    ggplot(aes(x = entry, y = corr, fill = entry)) + 
    facet_grid(env ~ cv, scales = "free_y") + 
    geom_boxplot(outlier.size = .5, linewidth = .4) +
    theme_bw() + theme(legend.position = "top", axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       strip.text.y = element_text(size = 8)) +
    labs(x = "Strategy", y = "Predictive ability", fill = "")  +
    scale_fill_manual(values = entrycol, labels = entrylab) +
    guides(fill = guide_legend(nrow = 1)),
  
  cvres.inter |>
    mutate(entry = factor(entry, levels = c("BLUE", "dBLUP", "dABLUP", "dABLUP_adj")))|> #reframe(corr = mean(corr), .by = c("entry",'pool','cv','repl')) |> 
    ggplot(aes(x = entry, y = mspe, fill = entry)) + 
    facet_grid(env ~ cv, scales = "free_y") + 
    geom_boxplot(outlier.size = .5, linewidth = .4) +
    theme_bw() + theme(legend.position = "top", axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       strip.text.y = element_text(size = 8)) +
    labs(x = "Strategy", y = "MSPE", fill = "") +
    scale_fill_manual(values = entrycol, labels = entrylab) +
    guides(fill = guide_legend(nrow = 1)),
  common.legend = TRUE, labels = "AUTO", nrow = 1
)


# True validation ---------------------------------------------------------
rm(list=ls())
load(file = "Saves/rval_pred.RDA")
load(file = "Saves/ptmods.RDA")

rval_across = do.call(rbind, lapply(blup.rval, function(x){
  left_join(
    x |> reframe(ots = mean(solution), .by = c('geno')),
    ptmods$blup |> reframe(pts = mean(solution), .by = 'geno')
  )
})) |> rownames_to_column('entry') |> mutate_at('entry', str_replace, '\\..*','') |> 
  mutate(entry = factor(
    case_when(entry == "dABLUP_adj" ~ "2S-dABLUP*",
              entry == "BLUE" ~ "2S-BLUE",
              entry == 'dBLUP' ~ "2S-dBLUP",
              entry == "dABLUP" ~ '2S-dABLUP'),
    levels = c("2S-BLUE", "2S-dBLUP", "2S-dABLUP", "2S-dABLUP*")
  ))

col_lab = c("Predictive ability", "MSPE")
names(col_lab) = c('corr', 'mspe')

ggarrange(
  do.call(rbind, lapply(split(rval_across, rval_across$entry), function(x) 
    data.frame(corr = cor(x$ots, x$pts, use = "na.or.complete"),
               mspe = mean((x$pts - x$ots)^2, na.rm = TRUE),
               entry = unique(x$entry)))) |> 
    pivot_longer(corr:mspe) |>
    ggplot(aes(x = entry, y = value)) + 
    facet_wrap(name~., scales = 'free_y', labeller = labeller(.cols = col_lab)) +
    geom_segment(aes(x = entry, xend = entry,
                     y = 0, yend = value), linewidth = 1) +
    geom_point(size = 3, color = 'black') + 
    theme_bw() + theme(axis.title = element_blank()),
  rval_across |>
    ggplot(aes(x = ots, y = pts)) + 
    facet_wrap(.~entry, scales = 'free_x', nrow = 1) +
    geom_point(alpha = .7) + 
    labs(x = "Predicted GEBVs using OTs", y = "Actual GEBVs from PYTs") + 
    theme_bw() +
    # theme(axis.title = element_blank()) +
    geom_smooth(method = 'lm', se = FALSE, colour = 'tomato'),
  nrow = 2, labels = 'AUTO'
)

# Training set optimization -----------------------------------------------
source("Codes/Preparations.R")

snp = snp[which(rownames(snp) %in% rownames(Gmat)),]
snp = snp[-duplicated(rownames(snp)),]
snpA = snp[which(rownames(snp) %in% unique(pheno[which(pheno$pool == "A"), c('geno')])),]
snpB = snp[which(rownames(snp) %in% unique(pheno[which(pheno$pool == "B"), c('geno')])),]

snp_scale = scale(snp, center = TRUE, scale = FALSE)
snpA_scale = scale(snpA, center = TRUE, scale = FALSE)
snpB_scale = scale(snpB, center = TRUE, scale = FALSE)

MlM = svd(snpA_scale)
L_A = snpA_scale %*% MlM$v; rm(snpA_scale)
MlM = svd(snpB_scale)
L_B = snpB_scale %*% MlM$v; rm(MlM, snpB_scale)

trsize_A = ceiling(nrow(L_A) * seq(.5, .8, .05))
trsize_B = ceiling(nrow(L_B) * seq(.5, .8, .05))

load("Saves/trt_opt_poolA.RDA")
load("Saves/trt_opt_poolB.RDA")

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

freqdfA = na.exclude(data.frame(geno = rownames(L_A)) |>
  full_join(
    do.call(rbind, lapply(freqlistA, function(x) do.call(rbind, x))) |> 
      rownames_to_column('trt_size') |>
      mutate_at("trt_size", str_replace, '\\..*', ''),
    by = c('geno' = "Var1")
  ))

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

freqdfB = na.exclude(data.frame(geno = rownames(L_B)) |>
  full_join(
    do.call(rbind, lapply(freqlistB, function(x) do.call(rbind, x))) |> 
      rownames_to_column('trt_size') |>
      mutate_at("trt_size", str_replace, '\\..*', ''),
    by = c('geno' = "Var1")
  ))

trtlabsA = paste0(paste0(seq(.5, .8, .05) * 100, '% '), "(",
                  gsub("size_", '', unique(freqdfA$trt_size)), ')')
names(trtlabsA) = unique(freqdfA$trt_size)
trtlabsB = paste0(paste0(seq(.5, .8, .05) * 100, '% '), "(", 
                  gsub("size_", '', unique(freqdfB$trt_size)), ')')
names(trtlabsB) = unique(freqdfB$trt_size)
critlabs = c("cd" = "CD", 'pev' = "PEV", 'rscore' = "r-Score")
cricols = c('cd'='#e41a1c', 'pev' = '#377eb8', 'rscore' = '#4daf4a')

## Selection frequency in the optimized training set -----------------------
ggarrange(
  ggplot(data = freqdfA, aes(x = geno, y = Freq)) + 
    theme_bw() + theme(legend.position = 'top', axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()) +
    facet_grid(Var2 ~ trt_size, labeller = labeller(.cols = trtlabsA, .rows = critlabs)) +
    geom_col(aes(fill = factor(trtopt))) +
    ylim(0, 40) + 
    scale_fill_manual(values = c("1" = "purple4", "0" = 'darkgrey'), 
                      labels = c("1" = 'Selected', '0' = "Non-selected")) + 
    labs(x = "Genotypes", y = "Frequency in the optimized training set", 
         fill = "Final optimized training set"),
  ggplot(data = freqdfB, aes(x = geno, y = Freq)) + 
    theme_bw() + theme(legend.position = 'top', axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()) +
    facet_grid(Var2 ~ trt_size, labeller = labeller(.cols = trtlabsB, .rows = critlabs)) +
    geom_col(aes(fill = factor(trtopt))) +
    ylim(0, 40) + 
    scale_fill_manual(values = c("1" = "purple4", "0" = 'darkgrey'), 
                      labels = c("1" = 'Selected', '0' = "Non-selected")) + 
    labs(x = "Genotypes", y = "Frequency in the optimized training set", 
         fill = "Final optimized training set"),
  common.legend = TRUE, labels = "AUTO", nrow = 2
)

## CV ---------------------------------------------------------------------
load(file = 'Saves/opt_cv_A.RDA')
load(file = 'Saves/random_cv_A.RDA')

opt_cv_resA = do.call(rbind, lapply(opt_cv, function(x) do.call(rbind, lapply(x, function(x){
  data.frame(
    corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
    mspe = mean((x$ybar - x$gebv)^2, na.rm = TRUE)
  )
})) %>% rownames_to_column("crit"))) %>% rownames_to_column('trt_size') %>% 
  mutate_at('trt_size',str_replace, '\\..*','') |> 
  mutate_at("trt_size", str_replace, 'size_','')
optA = rbind(
  opt_cv_resA |> mutate(size = as.numeric(trt_size)) |> select(-trt_size),
  do.call(rbind, random_cv) |> mutate(crit = 'random')
) 

trtlabsA = paste0(paste0(seq(.5, .8, .05) * 100, '% '), "(", unique(opt_cv_resA$trt_size), ')')
names(trtlabsA) = unique(opt_cv_resA$trt_size)

load(file = 'Saves/opt_cv_B.RDA')
load(file = 'Saves/random_cv_B.RDA')

opt_cv_resB = do.call(rbind, lapply(opt_cv, function(x) do.call(rbind, lapply(x, function(x){
  data.frame(
    corr = cor(x$gebv, x$ybar, use = 'na.or.complete'),
    mspe = mean((x$ybar - x$gebv)^2, na.rm = TRUE)
  )
})) %>% rownames_to_column("crit"))) %>% rownames_to_column('trt_size') %>% 
  mutate_at('trt_size',str_replace, '\\..*','') |> 
  mutate_at("trt_size", str_replace, 'size_','')
optB = rbind(
  opt_cv_resB |> mutate(size = as.numeric(trt_size)) |> select(-trt_size),
  do.call(rbind, random_cv) |> mutate(crit = 'random')
) 
trtlabsB = paste0(paste0(seq(.5, .8, .05) * 100, '% '), "(", unique(opt_cv_resB$trt_size), ')')
names(trtlabsB) = unique(opt_cv_resB$trt_size)

critlabs = c("cd" = "CD", 'pev' = "PEV", 'rscore' = "r-Score")
cricols = c('cd'='#e41a1c', 'pev' = '#377eb8', 'rscore' = '#4daf4a')
opt = rbind(optB |> mutate(pool = "B"), optA |> mutate(pool = "A"))

ggplot() + 
  geom_point(data = subset(opt, subset = crit == "random"), 
             aes(x = factor(size), y = corr), position = position_jitter(),
             alpha = .3) + 
  geom_boxplot(data = subset(opt, subset = crit == "random"), 
               aes(x = factor(size), y = corr), alpha = .5) + 
  geom_point(data = subset(opt, subset = crit != "random"), 
             aes(x = factor(size), y = corr, shape = crit, color = crit), 
             size = 3, alpha = .9) + 
  scale_colour_manual(values = cricols, labels = critlabs) + 
  scale_shape_manual(labels = critlabs, values = c(15, 17, 19))+ 
  theme_bw()  + theme(legend.position = 'top', legend.title = element_blank(),
                      axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)) +
  labs(x = "Training set sizes", y = "Predictive ability") + 
  scale_x_discrete(labels = c(trtlabsB, trtlabsA)) + 
  facet_wrap(.~pool, scales = 'free_x', labeller = labeller(.cols = c("A" = "Pool A", "B" = "Pool B")))


ggarrange(
  optA |> reframe(corr=  mean(corr), .by = c(size, crit)) |> 
    mutate(crit = factor(crit, levels = c("cd", 'pev', 'rscore', 'random'))) |> 
    ggplot() + 
    geom_point(aes(x = factor(size), y = corr, shape = crit, color = crit), 
               size = 3, alpha = .8) + 
    geom_line(aes(x = factor(size), y = corr, group = crit, color = crit)) + 
    scale_colour_manual(values = c(cricols, random = "black"), 
                        labels = c(critlabs, "random" = 'Random')) + 
    scale_shape_manual(labels = c(critlabs, "random" = 'Random'),
                       values = c("cd" = 15, "pev" = 17, "rscore" = 19, 'random' = 18))+ 
    theme_bw()  + theme(legend.position = 'top', legend.title = element_blank()) +
    labs(x = "Training set sizes", y = "Predictive ability") + 
    scale_x_discrete(labels = trtlabsA),
  optB |> reframe(corr=  mean(corr), .by = c(size, crit)) |> 
    mutate(crit = factor(crit, levels = c("cd", 'pev', 'rscore', 'random'))) |> 
    ggplot() + 
    geom_point(aes(x = factor(size), y = corr, shape = crit, color = crit), 
               size = 3, alpha = .8) + 
    geom_line(aes(x = factor(size), y = corr, group = crit, color = crit)) + 
    scale_colour_manual(values = c(cricols, random = "black"), 
                        labels = c(critlabs, "random" = 'Random')) + 
    scale_shape_manual(labels = c(critlabs, "random" = 'Random'),
                       values = c("cd" = 15, "pev" = 17, "rscore" = 19, 'random' = 18))+ 
    theme_bw()  + theme(legend.position = 'top', legend.title = element_blank()) +
    labs(x = "Training set sizes", y = "Predictive ability") + 
    scale_x_discrete(labels = trtlabsB),
  labels = "AUTO", common.legend = TRUE, nrow = 2
)
