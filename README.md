# 🍠 Towards a sweetpotato genomic-enabled breeding: optimizing two-stage analysis of multi-environment augmented trials

This repository has the data and the analytical pipeline adopted by Chaves et al. (2026) (manuscript currently under review).

Although single-stage models are generally more efficient, there are contexts where they are difficult to fit, making two-stage models the most practical alternative. In this study, we used observational trials (OTs, implemented using augmented row–column designs) from sweetpotato breeding to evaluate different two-stage strategies, using the single-stage (SS) model as a benchmark. Given the lack of replication in OTs, we hypothesized that deregressed best linear unbiased predictions (dBLUPs) or pedigree-based dBLUPs (dABLUPs) would work more appropriately as inputs for second-stage models than best linear unbiased estimates (BLUEs). These comparisons were conducted within weighted models using either a diagonal weight matrix or the full weight matrix. Our analyses focused on storage root yield across six OTs in Uganda, derived from intra-genepool crosses (1,138 test clones in total). We compared the selection and prediction performance of genomic BLUP models under different two-stage strategies against the SS benchmark, and we also tested whether pool-specific genomic prediction models offered advantages over models trained with the complete dataset. For selection, differences among second-stage models were minor, with a slight advantage for those using dABLUPs as entries, combined with the full weight matrix. For prediction, however, the choice of weighting scheme had a greater impact on performance than the choice of entry. Using the complete dataset, differences between entries were marginal, but for pool-specific predictions, dABLUPs provided the best performance. Overall, if adopting a two-stage strategy for the analysis of augmented trials, we recommend using dABLUPs together with the full weight matrix.

## 🗝️ Key message

Using the full weight matrix and deregressed pedigree-based best linear unbiased predictions in second-stage models lead to selections and genomic predictions closer to those obtained using a single-stage model.

<img width="3160" height="2784" alt="gebvs_stg2ss" src="https://github.com/user-attachments/assets/223583e1-63fb-4557-adb6-d6e543b4a1a3" />

<img width="3160" height="2784" alt="CV_both" src="https://github.com/user-attachments/assets/acc4e1e6-9b9f-4679-9d57-1a5218ed4f82" />


## 🛠️ Repository structure

```
🌳 schaves_sweetpotato_two_stage/
├── 📂 Data/
│   ├── Phenotype.csv                   # Raw phenotypic data from multi-environment trials
│   ├── SNPs.RDS                        # Hexaploid SNP dosage matrix (0 to 6)
│   └── Pedigree.csv                    # Pedigree file
├── 📂 Codes/
│   ├── 1_preparations.R                # Raw data management
│   ├── 2_first_stage.R                 # Fitting single-environment models
│   ├── 3_second_stage.R                # Fitting second-stage multi-environment models
│   ├── 4_single_stage.R                # Fitting single-stage multi-environment models
│   ├── 5_cross_validations             # Cross-validations
│   └── 📂 Functions/        
│       ├── expvar.R                    # Calculates the percentage of explained variance in factor analytic mixed models
│       ├── fa_summa_rr.R               # Obtains rotated loadings and scores and builds the between-environment genetic variance-covariance and correlation matrices 
│       ├── measure_variances.R         # Partitions the genotype-by-environment interaction into different pars (by Daniel Tolhurst and Christian Werner)
│       └── upmod.R                     # Updates a model fitted using asreml until complete convergence
├── README.md                           # Project documentation
├── LICENSE                             # Project license
├── schaves_sweetpotato_two_stage.Rproj # R project file
└── .gitignore                          # Files excluded from version control
```

## 👨‍🏫 Citation

Soon.

## 📧 Contact

Please, feel free to contact for any enquiries, suggestions or issues:

*saulochaves@usp.br*  
*g.pereira@ufv.br*


