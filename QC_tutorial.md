# Standard genotyping data QC

The standard QC procedure starts with samples QC and then variants QC. Analyses made in this document use plink and R and the examples assume you have genetic data in plink binary format bed/bim/fam.

Essential steps are described below

- [Standard genotyping data QC](#standard-genotyping-data-qc)
  - [Samples QC](#samples-qc)
    - [1. Missingness and het rate](#1-missingness-and-het-rate)
      - [Compute missing calls per sample](#compute-missing-calls-per-sample)
      - [Calculate heterozygosity](#calculate-heterozygosity)
      - [Plot and decide filtering thresholds](#plot-and-decide-filtering-thresholds)
    - [2. Samples IBD](#2-samples-ibd)
      - [Prepare a subset of indep SNPs](#prepare-a-subset-of-indep-snps)
      - [Compute kinship and IBD](#compute-kinship-and-ibd)
      - [Find highly related individuals](#find-highly-related-individuals)
    - [3. Sex check](#3-sex-check)
    - [4. Remove samples failing QC](#4-remove-samples-failing-qc)
  - [Variants QC](#variants-qc)
    - [1. Separate chrY](#1-separate-chry)
    - [2. Compute missing rate, allele frequency and HWE test](#2-compute-missing-rate-allele-frequency-and-hwe-test)
    - [(Optional) Differential missing rate between cases and controls](#optional-differential-missing-rate-between-cases-and-controls)
    - [3. Generate plots](#3-generate-plots)
    - [4. Apply filters for SNPs](#4-apply-filters-for-snps)
  - [Conclusion](#conclusion)

## Samples QC

### 1. Missingness and het rate

These are usually done excluding the sex chromosomes

#### Compute missing calls per sample

```bash
plink --bfile input.file --missing --chr 1-22 --allow-no-sex --out QC_step1.miss
```

#### Calculate heterozygosity

```bash
plink --bfile input.file --het --chr 1-22 --allow-no-sex --out QC_step1.het
```

#### Plot and decide filtering thresholds

Now we need to plot the metrics computed above and estimate where to put our threasholds. This can be done easily in R starting from the files generated above.

In general, too high missingness can indicate a poor quality sample, while high het rate usually indicates sample contamination. The exact distribution of het rate depends on cohort ancestry and low het rate is usually associate with high level of autozygosity (inbreeding).

A good starting point for filtering is:

- per sample missing 0.05
- het rate betwen 0.1 and 0.4

But it is always better to adjust them after plot inspection.

```R
#load data
het <- read.table("QC_step1.het.het", header=T)
mis <- read.table("QC_step1.miss.imiss", header=T)

#calculate heterozigosity rate
mishet <- data.frame(FID=het$FID, IID=het$IID, het.rate=(het$N.NM - het$O.HOM)/het$N.NM, mis.rate=mis$F_MISS)

#plot
plot(y=mishet$het.rate, x=mishet$mis.rate, ylab="Heterozygosity rate", xlab="Proportion of missing genotypes")

#Examine plot and decide appropriate threshold for filtering het and miss
miss.threshold = 0.05
het.nsd.threshold = 3 #N std from het rate mean for filtering

het.sd <- sd(mishet$het.rate, na.rm=T)
het.mean <- mean(mishet$het.rate, na.rm=T)
lower.het.threshold <- het.mean - (3 * het.sd)
upper.het.threshold <- het.mean + (3 * het.sd)


#Save plot image with desired thresholds 
png("QC_step1_mis_het.png", res=1200, height=5, width=5, units="in")
plot(y=mishet$het.rate, x=mishet$mis.rate, ylab="Heterozygosity rate", xlab="Proportion of missing genotypes")
abline(v=miss.threshold, lty=2) #vertical dashed line
abline(h=lower.het.threshold, lty=2) #horizontal dashed line
abline(h=upper.het.threshold, lty=2) #horizontal dashed line
dev.off()

#save filtered individuals to file
write.table(mishet[mishet$mis.rate > miss.threshold,], file="fail_mis_qc.txt", row.names=F, quote=F)
write.table(mishet[mishet$het.rate < lower.het.threshold | mishet$het.rate > upper.het.threshold,], file="fail_het_qc.txt", row.names=F, quote=F)
```

### 2. Samples IBD

Here we aim to remove individuals with unexpected relationships.]

#### Prepare a subset of indep SNPs

To speed up IBD computation we prepare a small subset of filtered common SNPs and then prune them.

```bash
#create a subset of frequent SNPs, with low missing rate
plink --bfile input.file --chr 1-22 --maf 0.25 --geno 0.05 --hwe 1e-6 --allow-no-sex --make-bed --out QC_step2_frequent

#perform pruning
plink --bfile QC_step2_frequent --indep-pairwise 50 5 0.2 --allow-no-sex --out QC_step2_prunedSNPs

#generate IBD report from pruned SNPs
plink --bfile QC_step2_frequent --extract QC_step2_prunedSNPs.prune.in --genome --allow-no-sex --out QC_step2_pruned
```

#### Compute kinship and IBD

Relationships between individuals can be estimated using plink `--genome` or king `--related` options. The latter is usually more efficient for large cohorts.

```bash
#Using plink
plink --bfile QC_step2_pruned \
  --genome \
  --out QC_step2.ibd

#Using king
king -b QC_step2_pruned.bed \
  --related --cpus 8 \
  --prefix QC_step2.ibd
```

Plink computes PI_HAT as a metric for relatedness. This can be interpreted as the fraction of shared genome between 2 individuals, thus the expected value is 1, 0.5, 0.25, 0.125 for duplicate/MZ twin, 1st-degree, 2nd-degree and 3rd-degree relationships, respectively.

King estimates a kinship coefficient that range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] for duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively. Note that only pairs with kingship > 0.177 are emitted by default when using `--related` as described above.

#### Find highly related individuals

Then we can use R to compute IBD, plot and filter individuals based on the chosen threshold. Usually, we want to remove either duplicated samples or close related individuals (1st degree). One can adjust thresholds value according to the need the specific project.

In case of plink genome data, you can use something like this. Here we will remove only samples with very high PI_HAT > 0.5 that are likely to be duplicated or 1st degree relatives.

```R
#read file
genome <- read.table("QC_step2.ibd.genome", header=T, as.is=T)

#select pairs of individuals with PI_HAT > 0.8 (likely duplicated samples)
pihat_max <- 0.5
genome <- genome[genome$PI_HAT > pihat_max,]

#compute mean(IBD), var(IBD), SE(IBD)
mean.IBD <- 0*genome$Z0 + 1*genome$Z1 + 2*genome$Z2
var.IBD <- ((0-mean.IBD)^2)*genome$Z0 + ((1-mean.IBD)^2)*genome$Z1 + ((2-mean.IBD)^2)*genome$Z2
se.IBD <- sqrt(var.IBD)

#save IBD plot for individuals with high PI_HAT
png("QC_step2_IBD.png", res=1200, height=5, width=5, units="in")
plot(mean.IBD, se.IBD, xlab="Mean IBD", ylab="SE IBD")
dev.off()

#Find individuals with mean IBD > 1 and for each couple filter out the sample with higher missing genotypes
high_ibd <- genome[mean.IBD > 1, ]
mis <- read.table("QC_step1.miss.imiss", header=T)
high_ibd <- merge(high_ibd, mis[,c(1,2,6)], by.x=c("FID1","IID1"), by.y=c("FID","IID")) #Suppose unique ids are in FID field
colnames(high_ibd)[ncol(high_ibd)] <- "F_MISS1"
high_ibd <- merge(high_ibd, mis[,c(1,2,6)], by.x=c("FID2","IID2"), by.y=c("FID","IID"))
colnames(high_ibd)[ncol(high_ibd)] <- "F_MISS2"
fail_ibd_qc <- data.frame(FID=character(), IID=character(), stringsAsFactors=F)

for (n in 1:nrow(high_ibd)) {
  if (high_ibd$F_MISS1[n] > high_ibd$F_MISS2[n]) {
    fail_ibd_qc[nrow(fail_ibd_qc)+1,] <- c(high_ibd$FID1[n], high_ibd$IID1[n])
  } else {
    fail_ibd_qc[nrow(fail_ibd_qc)+1,] <- c(high_ibd$FID2[n], high_ibd$IID2[n])
  }
}

#Save list of individuals failing IBD QC
write.table(fail_ibd_qc, "fail_ibd_qc.txt", row.names=F, quote=F)
```

In case of king related files you can use this. The default output using `--related` contains only closely related samples with kingship > 0.177 (1st degree)

```R
kinship <- read.table("QC_step2.ibd.kin0", header=T, as.is=T)
ggplot(kinship, aes(x=Kinship)) + geom_histogram(bin=100) + theme_bw()
ggplot(kinship, aes(x=InfType)) + geom_bar() + scale_y_sqrt() + theme_bw()

mis <- read.table("QC_step1.miss.imiss", header=T)
high_ibd <- merge(kinship, mis[,c(1,2,6)], by.x=c("FID1","ID1"), by.y=c("FID","IID")) #Suppose unique ids are in FID field
colnames(high_ibd)[ncol(high_ibd)] <- "F_MISS1"
high_ibd <- merge(high_ibd, mis[,c(1,2,6)], by.x=c("FID2","ID2"), by.y=c("FID","IID"))
colnames(high_ibd)[ncol(high_ibd)] <- "F_MISS2"
fail_ibd_qc <- data.frame(FID=character(), IID=character(), stringsAsFactors=F)

for (n in 1:nrow(high_ibd)) {
  if (high_ibd$F_MISS1[n] > high_ibd$F_MISS2[n]) {
    fail_ibd_qc[nrow(fail_ibd_qc)+1,] <- c(high_ibd$FID1[n], high_ibd$IID1[n])
  } else {
    fail_ibd_qc[nrow(fail_ibd_qc)+1,] <- c(high_ibd$FID2[n], high_ibd$IID2[n])
  }
}

#Save list of individuals failing IBD QC
write.table(fail_ibd_qc, "fail_ibd_qc.txt", row.names=F, quote=F)
```

### 3. Sex check

To perform sex check using plink we first need to ensure that pseudo-autosomal regions have been removed. We can use `--split-x` for this. Possible genome builds are `hg18`, `hg19` and `hg38`

```bash
plink --bfile input.file \
  --split-x hg38 no-fail \
  --make-bed --out input.file.XY
```

Then use plink to perform the sex check. We suggest to use a mild threshold (0.5, 0.8).

```bash
plink --bfile input.file.XY --sex-check 0.5 0.8 --out sex-check
```

In the resulting sex check file, status is 'OK' if PEDSEX and SNPSEX match and are nonzero, 'PROBLEM' otherwise.
Keep in mind that this is based on the F stats and you may need to adjust thresholds. Usually it is better to run `--check-sex` once, eyeballing the distribution of F estimates. There should be a clear gap between a very tight male clump at the right side of the distribution and the females everywhere else. After looking at the distribution one can rerun with parameters corresponding to the empirical gap. For example, same samples can end up as undetermined based on SNP data and inspection of actual F values may help claryfing their status.

### 4. Remove samples failing QC

We can now remove the samples failing any of the above QC steps

```bash
#merge single fail_qc lists
(head -n1 fail_het_qc.txt | cut -d" " -f1,2 & cut -d" " -f1,2 fail_* | grep -v "FID IID" | sort -u) > fail_qc.txt

#remove samples failing qc
plink --bfile input.file --remove fail_qc.txt --allow-no-sex --make-bed --out QCed_samples
```

After that we end up with a new dataset that we call `QCed_samples`. Using this we can proceed to variant level QC.

## Variants QC

### 1. Separate chrY

Variants on chrY will be present only in males and they will have a high missing rate if we consider them across the whole cohort. Thus we will separate variants on chrY to a distinct file and apply QC to variants on autosomes + chrX.

```bash
plink --bfile QCed_samples --chr Y --out QCed_dataset.chrY
plink --bfile QCed_samples --not-chr Y --out QCed_samples.chr1-23MT
```

### 2. Compute missing rate, allele frequency and HWE test

Starting from a sample QCed dataset we first compute missing rate for each variant and the allele frequency

```bash
#compute SNPs missing rate
plink --bfile QCed_samples.chr1-23 --missing --allow-no-sex --out QC_step3.mis

#compute allele freq
plink --bfile QCed_samples.chr1-23 --freq --allow-no-sex --out QC_step3.freq

# compute HWE and HWE violation test
plink --bfile QCed_samples.chr1-23 --hardy --allow-no-sex --out QC_step3.hwe
```

### (Optional) Differential missing rate between cases and controls

If we are working with a case-control cohort, it is suggested to compare the difference in missing rate between cases and control to remove SNPs with largely different missing rate between the 2 groups. plink can compute a test for differential missingness using `--test-missing`

```bash
#suppose phenos are encoded in your files, otherwise provide pheno with --pheno option
plink --bfile QCed_samples.chr1-23 --test-missing --allow-no-sex --out QC_step3.diffmiss
```

### 3. Generate plots

Use R to generate variants QC plots and eventually list of SNPs failing the differential missingness test.

```R
#create list of SNPs filtered for differential call rate in cases/controls
diffmiss <- read.table("QC_step3.diffmiss.missing", header=T, as.is=T)
diffmiss <- diffmiss[diffmiss$P < 1e-5,]
write.table(diffmiss$SNP, "fail_SNPs_diffmiss.txt", row.names=F, col.names=F, quote=F)

#create plots for MAF and missing rate
lmis <- read.table("QC_step3.mis.lmiss", header=T)
png("SNPs_lmiss.png", res=1200, width=5, height=5, units="in")
plot(log10(lmis$F_MISS), ylab="Number of SNPs", xlab="Fraction missing genotypes (log10)", main="Fraction of missing data")
abline(v=log10(0.05), lty=2)
dev.off()

freq <- read.table("QC_step3.freq.frq", header=T)
png("SNPs_MAF.png", res=1200, width=5, height=5, units="in")
hist(freq$MAF, ylab="Number of SNPs", xlab="MAF", main="Minor allele frequencies")
abline(v=0.01, lty=2)
dev.off()

hwe <- read.table("QC_step3.hwe.hwe", header=T)
png("SNPs_HWE.png", res=1200, width=5, height=5, units="in")
hist(round(-log10(hwe$HWE), 4), ylab="Number of SNPs", xlab="-log(HWE test P)", main="HWE test")
abline(v=15, lty=2)
dev.off()
```

### 4. Apply filters for SNPs

Based on inspection of the plots generated above we can finally set filters for SNPs. Usually, we can use plink to remove SNPs with missing genotypes > 5% (`--geno`) and severely violating HWE (`--hwe`). If we identified SNPs with differential missingness across cases and controls we can remove them as well (`--exclude`). An example command is reporte below

```bash
#Filter out unwanted SNPs
plink --bfile QCed_samples.chr1-23MT \
    --hwe 1e-15 --geno 0.05 \
    --exclude fail_SNPs_diffmiss.txt \ #Optional
    --allow-no-sex --make-bed --out QCed_dataset.chr1-23MT
```

## Conclusion

At the end of the QC procedure, we end up with a `QCed_dataset.chr1-23MT` in bed/bim/fam format that contains autosomes, the X chromosome and MT variants ready for downstream analyses.

Additionally, variants from the Y chromosome are stored separately in `QCed_dataset.chrY`.
