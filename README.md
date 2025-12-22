# Tractor-Mix


Tractor-Mix is a GWAS method tailored for admixed populations with related individuals. It shares the core intuition of [Tractor](https://github.com/Atkinson-Lab/Tractor) by assigning genetic markers based on local ancestry and jointly modeling them in a regression framework. To address sample relatedness, we developed a linear/logistic mixed model incorporating an estimated GRM, allowing us to generate joint and ancestry-specific p-values, as well as effect size estimates.

Tractor-Mix consists of 6 steps, namely:  
* Phasing and local ancestry inference  
* Deconvolute genetic markers by local ancestry  
* Estimate covariates with PC-Air (or basic PC)
* Estimate GRM/kinship with PC-Relate  (or basic GRM)
* Fit the null model with [GMMAT](https://github.com/hanchenphd/GMMAT)  
* Run Tractor-Mix

Note: The current implementation can handle a cohort with ~15000 individuals; a larger cohort might take days to finish.

![Image Title](pipeline.png)


&nbsp;  
&nbsp;  


### Dependencies: 

R > 4.0 with `GMMAT`, `Matrix`, `data.table`, `doParallel`, `foreach`, `dplyr`, `abind` pre-installed


### Step1: Phasing and local ancestry inference  

Similar to Tractor, Tractor-Mix relies on the effective performance of phasing and local ancestry inference. We have developed a tutorial on using Shapeit2/Shapeit5 and Rfmix2; for more details, please refer to [this page](https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Rfmix.md).

When selecting a reference panel, it’s important to choose continental or subcontinental populations that closely match the query admixture. For example, in an African-American cohort, using all continental AFR and EUR populations from the 1000 Genomes Project would be advisable as a reference panel. Excessive ancestry calls can diminish statistical power. For instance, a cohort with [AFR, EUR, AMR] = [60%, 38%, 2%] is not ideal for Tractor-Mix due to the low 2% contribution from AMR. In such cases, treating the cohort as a 2-way admixture is more appropriate.


&nbsp;  
&nbsp;  

### Step2: Deconvolute genetic markers by local ancestry  

To allocate alleles to each local ancestry background, we will use the `extract_tracts.py` script from Tractor. For detailed instructions, please refer to [this page](https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Extract.md). Briefly, you can execute the following code for this step:

```
python3 extract_tracts.py \
  --vcf file.phased.vcf.gz \
  --msp file.deconvoluted.msp.tsv \
  --num-ancs 2
```
For a 2-way admixture, this process will generate six files in total (2 VCFs, 2 hapcount.txt, and 2 dosage.txt files). In our Tractor-Mix model, we will only use the 2 dosage.txt files.


&nbsp;  
&nbsp;  

### Step3: Estimate covariates with PC-Air [or simple PC]

In the linear/logistic mixed model, principal components (PCs) or admixture proportions are typically used to account for population stratification. For cohorts with related samples, a common approach is to first partition the samples into independent and related subsets. We construct the PC space using the independent samples and then project the related individuals into this space, ensuring that the principal components are not influenced by family structure. This method is implemented in the R package [GENESIS](https://github.com/UW-GAC/GENESIS), which can be installed via Bioconductor. A tutorial on running PC-Air can be found [here](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html).

Briefly, if you have VCF files, you might perform LD pruning using PLINK and then convert the `bim/bam/fam` files to GDS format with `SNPRelate::snpgdsBED2GDS`. In addition to genotype data, GENESIS requires a KING-robust estimate as input. You can compute the KING-robust matrix using `SNPRelate`, as instructed in the GENESIS vignette, or use PLINK2 with the `plink2 --make-king square` argument.

Alternatively, you may use PLINK to calculate the PC.

**[Note: In our simulations and the benchmark from the [GENESIS paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7904076/#sup1), the difference between raw PCA vs PC-Air is minor.]**



&nbsp;  
&nbsp;  

### Step4: Estimate GRM/kinship with PC-Relate  [or standard GRM]

After completing PC-Air, you can run PC-Relate to obtain the estimated kinship matrix, using several top PCs from PC-Air (refer to [the same tutorial](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html)). Please note that GENESIS may occasionally scramble the order of samples in the output file, especially if your sample IDs are numerical, so it's important to ensure the order is correct.

To improve computational efficiency, we recommend creating a sparse GRM using the following code:
```
library(Matrix)
PCRelatemat = pcrelateToMatrix(pcrelate_res, sample.include = iids, scaleKin = 2)[iids,iids]
PCRelatemat_sparse = PCRelatemat
# mask values < 0.05
PCRelatemat_sparse[PCRelatemat_sparse < 0.05] = 0
PCRelatemat_sparse = as(PCRelatemat_sparse, "sparseMatrix") 
```

Alternatively, you may use PLINK/GCTA/SAIGE to calculate the GRM. But you may need to reformat the GRM to be a square matrix, so that it can be read by GMMAT. You shall consider masking a GRM, and convert it to `sparseMatrix`. This procedure often largely enhance the computational efficiency. 

**[Note: In our simulations and the benchmark from the [GENESIS paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7904076/#sup1), the difference between raw GRM vs PC-Relate is minor]**

&nbsp;  
&nbsp;  

### Step5: Fit the null model with [GMMAT](https://github.com/hanchenphd/GMMAT)  


GMMAT can be installed using `install.packages("GMMAT")`. To fit the null model, use the following commands:

For a continuous outcome:
```
# continuous
Model_Null =   glmmkin(fixed = Pheno ~  Covariates, 
                       data = YourDf, id = "ID", kins = GRM, 
                       family = gaussian())

# dichotomous
Model_Null = glmmkin(fixed = Pheno ~ Covariates, 
                     data = YourDf, id = "ID", kins = GRM, 
                     family = binomial())
```

In essence, we fit a linear or logistic mixed model that includes covariates as the fixed effect term and the GRM as the random effect term. This null model will later be used for the genome-wide scan.

*Please ensure that the order of IDs in your covariate file matches the order in the GRM and dosage files.*


&nbsp;  
&nbsp;  

### Step6: Run Tractor-Mix (unconditional or conditional)

Now you should have the null model and ancestry-specific genotype `dosage` files. You can run unconditional Tractor-Mix with `TractorMix.score`. You can also run conditional analysis by including the `hapcount` file in `TractorMix.score_cond`.
```
# unconditional
source("TractorMix.score.R")
TractorMix.score(obj = Model_Null, 
                 infiles = c("Admix.anc0.dosage.txt", "Admix.anc1.dosage.txt"),
                 outfiles = "result.tsv", 
                 AC_threshold = 50)
 
# conditional      
source("TractorMix.score_cond.R")
TractorMix.score_cond(obj = Model_Null, 
                 infiles_geno = c("Admix.anc0.dosage.txt", "Admix.anc1.dosage.txt"),
                 infiles_la = c("Admix.anc0.hapcount.txt"),
                 outfiles = "result.tsv", 
                 AC_threshold = 50)
```


From TractorMix, you should obtain:  
* `P`: a joint p-value
* `Eff_anc0`, `Eff_anc1`: ancestry-specific effect size estimates
* `SE_anc0`, `SE_anc1`: ancestry-specific standard error
* `Pval_anc0`, `Pval_anc1`: ancestry-specific p values
* `include_anc0`, `include_anc1`: indicate if ancestry-specific dosages are included (to prevent false positive inflation due to low MAC)



In our UKBB analysis with 9,000 samples and 8.5M variants, it took approximately 3 hours per trait with 22 chromosomes run in parallel on a typical high performance computing cluster setup (using sparse GRM, 16GB RAM, 8 CPUs per chromosome). 


&nbsp;  
&nbsp;  

### Optional: Wald test of Tractor-Mix (TBD)

Wald test should only be applied to a smaller set of markers (e.g. less than 100 variants). 







