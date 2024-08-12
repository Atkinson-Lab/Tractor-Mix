# Tractor-Mix

TractorMix is a GWAS method designed for admixed populations with related individuals. TractorMix inherits similar intuition as [Tractor](https://github.com/Atkinson-Lab/Tractor), where we allocate genetic markers by local ancestry, and jointly model them in a regression setting. To account for sample relatedness, we built a linear/logistic mixed model with an estimated GRM, and produce joint/ancestry-specific p-values and effect size estimates.

TractorMix consists of 6 steps, namely:  
* Phasing and local ancestry inference  
* Deconvolute genetic markers by local ancestry  
* Estimate covariates with PC-Air  
* Estimate GRM/kinship with PC-Relate  
* Fit the null model with [GMMAT](https://github.com/hanchenphd/GMMAT)  
* Run Tractor-Mix


![Image Title](pipeline.png)


&nbsp;  
&nbsp;  

### Step1: Phasing and local ancestry inference  

Similar to Tractor, Tractor-Mix requires good performances of phasing and local ancestry inference. We have developed a tutorial on how to use Shapeit2/Shapet5 and Rfmix2, please refer to [this page](https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Rfmix.md) for more details.

For the reference panel, users may want to use continental/subcontinental populations that accurately match the query admixture. (e.g. for the African-American cohort, one may want to use all continental AFR and EUR populations from 1000 Genome projects as the reference panel). Excessive ancestry calls can reduce statistical power. For instance, a cohort with [AFR, EUR, AMR] = [60%, 38%, 2%] is not preferable for Tractor-Mix since AMR only composes 2%. Instead, calling such cohort 2-way admixture is more appropriate.


&nbsp;  
&nbsp;  

### Step2: Deconvolute genetic markers by local ancestry  

We will use `extract_tracts.py` code from Tractor to allocate alleles to each local ancestry background. Please refer to [this page](https://github.com/Atkinson-Lab/Tractor-tutorial/blob/main/Extract.md) for more details. Briefly, you can run the following code for this step:
```
python3 extract_tracts.py \
  --vcf file.phased.vcf.gz \
  --msp file.deconvoluted.msp.tsv \
  --num-ancs 2
```
For 2-way admixture, we will obtain 6 files in total (2 vcfs, 2 hapcount.txt, 2 dosage.txt). We will only use 2 dosage.txt files in our TractorMix model. 


&nbsp;  
&nbsp;  

### Step3: Estimate covariates with PC-Air  

In the linear/logistic mixed model, we typically use PC or admixture proportion to capture population stratification. For a cohort with related samples, a common practice is to first partition samples into independent subsets and relative subsets. We construct the PC space with independent samples, then project relatives into the space. Through this approach, principal components will not be driven by family structure. This approach has been implemented in the R package [GENESIS](https://github.com/UW-GAC/GENESIS), which can be installed through Bioconductor. The tutorial on running PC-Air can be found [here](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html). Briefly, if you have vcf files in hand, you may need to convert them into GDS file formats so that they can be processed by GENESIS (`vcf ---> PLINK ---> GDS`). GENESIS then performs LD prune, partition samples, and runs PC-Air. Once finished running PC-Air, one may want to save the results as a tsv file.

**[Note: In our current theoretical simulations and [GENESIS paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7904076/#sup1), the difference between raw PCA vs PC-Air is minor]**

&nbsp;  
&nbsp;  

### Step4: Estimate GRM/kinship with PC-Relate  

Once finish running PC-Air, once can run PC-Relate (refer to [the same tutorial](http://bioconductor.org/packages/release/bioc/vignettes/GENESIS/inst/doc/pcair.html)) to obtain the estimated kinship matrix. Please note that occasionally, GENESIS can scramble the order of samples in the output file, please make sure the order is correct. 

**[Note: In our current theoretical simulations and [GENESIS paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7904076/#sup1), the difference between raw GRM vs PC-Relate is minor]**

&nbsp;  
&nbsp;  

### Step5: Fit the null model with [GMMAT](https://github.com/hanchenphd/GMMAT)  

GMMAT can be installed with `install.packages("GMMAT")`. We will use the following command to fit the null model: 
```
# continuous:
Model_Null =   glmmkin(fixed = Pheno ~  Covariates, 
                       data = YourDf, id = "ID", kins = GRM, 
                       family = gaussian())

# dichotomous:
Model_Null =   glmmkin(fixed = Pheno ~  Covariates, 
                       data = YourDf, id = "ID", kins = GRM, 
                       family = binomial())
```

In essence, we fit a linear/logistic mixed model that takes covariates as the fixed effect term, and takes GRM as the random effect term. This is the first step of Rao's score test, where we fit the model under the null hypothesis and obtain its gradient and information matrix. 

For now, please ensure that the order of ID in your covariate file matches with the GRM order; Please also ensure there is no missing value in GRM, phenotype, or covariates. I will try to implement a better version of TractorMix to handle this down the road. 


&nbsp;  
&nbsp;  

### Step6: Run Tractor-Mix  

Now you should have the null model and ancestry-specific genotype dosage files. You can run TractorMix with:
```
source("TractorMix.score.R")

# continuous
TractorMix.score(obj = Model_Null, 
                 infiles = c("Genotype.anc0.hapcount.txt", "Genotype.anc1.hapcount.txt"),
                 outfiles = "result.tsv" )
                 
# dichotomous (need to specify a threshold to filter out variants with low ancestry-specific allele counts)
TractorMix.score(obj = Model_Null, 
                 infiles = c("Genotype.anc0.hapcount.txt", "Genotype.anc1.hapcount.txt"),
                 outfiles = "result.tsv", 
                 AC_threshold = 50)
```

From TractorMix, you should obtain:  
* a joint p-value (to reject the null hypothesis H0: beta(anc0) = beta(anc1) = 0)  
* ancestry-specific effect size estimates (Eff0, Eff1)
* ancestry-specific p values (PVal0, PVal1)

In our UKBB analysis with 9,000 samples and 8.5M variants, it took approximately 10 hours per trait with 22 chromosomes run in parallel on a typical high performance computing cluster setup (using sparse GRM, 32GB RAM, 8 CPUs per chromosomes). 




### Optional: Wald test of Tractor-Mix 

Wald test should only be applied to a smaller set of markers (e.g. less than 100 variants). 







