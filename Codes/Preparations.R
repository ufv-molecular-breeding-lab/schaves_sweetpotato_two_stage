rm(list=ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(asreml)
library(AGHmatrix)
library(desplot)

# Phenotypic data ---------------------------------------------------------
pheno = read.csv(file = "Data/Phenotype.csv")
head(pheno)
str(pheno)

## Overview ----------------------------------------------------------------
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


# Pedigree ----------------------------------------------------------------
ped = read.csv(file = "Data/Pedigree.csv")
head(ped)
str(ped)

Amat = Amatrix(data = ped[,-4], ploidy = 6, dominance = FALSE)

# Genomic data ------------------------------------------------------------
snp = readRDS(file = 'Data/SNPs.RDS')
table(snp)
genofreq = apply(snp, 2, function(x) table(x))
monom = names(which(do.call(c, lapply(genofreq, function(x) length(x) == 1))))
snp = snp[,-which(colnames(snp) %in% monom)]
genofreq = genofreq[which(!names(genofreq) %in% monom)]
dose = (0:6)
al = NULL
for (i in dose)
  al[i + 1] = paste(paste(rep("A", i), collapse = ""), paste(rep("a", (6 - i)), collapse = ""), sep = "")
names(dose) = al
maf = lapply(genofreq, function(x){
  A = sum(dose[which(dose %in% names(x))]/6 * x/sum(x))
  a = 1 - A
  min(a,A)
})
which(do.call(c, lapply(maf, function(x) x < 0.05)))

Gmat = Gmatrix(SNPmatrix = snp[-which(duplicated(rownames(snp))),],
               method = "VanRaden", ploidy = 6)
Gmat = Gmat[which(rownames(Gmat) %in% c(rownames(Amat))),
            which(colnames(Gmat) %in% c(colnames(Amat)))]

# Wrap-up -----------------------------------------------------------------
rm(list=ls()[which(!ls() %in% c("ped", "Gmat", "snp", "pheno", "Amat", "pedigree"))])

