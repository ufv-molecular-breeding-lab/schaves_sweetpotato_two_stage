rm(list=ls())

# Preps -----------------------------------------------------------------
source("Codes/Preparations.R")

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
PC = FactoMineR::PCA(X = Gmat, graph = FALSE, scale.unit = FALSE)
pc.plot = data.frame(
  geno = names(PC$ind$coord[,1]),
  pc1 = PC$ind$coord[,1],
  pc2 = PC$ind$coord[,2], row.names = NULL
) |> left_join(unique(pheno[,c("geno",'pool','status')]), by = "geno")

pc.plot |> filter(status != "Check") |>
  ggplot(aes(x = pc1, y = pc2)) +
  geom_vline(aes(xintercept = 0), linetype = "dashed") +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  geom_point(
    aes(color = pool, shape = status),
    size = 2,
    alpha = .8
  ) +
  ggrepel::geom_label_repel(
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
    y = paste0('PC2 (', round(PC$eig[2, 2], 2), '%)')
  ) +
  scale_color_manual("Pool", values = c('#1b9e77', '#d95f02')) +
  scale_fill_manual("Pool", values = c('#1b9e77', '#d95f02'))  +
  scale_shape_manual("", values = c(15, 17)) +
  theme_bw() + theme(legend.position = 'bottom')

