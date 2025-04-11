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

