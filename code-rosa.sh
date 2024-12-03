

# reformat phenos
cat gwas-phenotypes.csv | tr ',' ' ' > divorce_phen.txt

cp -r /tsd/p805/data/durable/users/torkildl/genetics-divorce/divorce_phen.txt /cluster/projects/p805/rosac/divorce/input/
 # edit phenofile (no header in gcta format; 
cd /cluster/projects/p805/rosac/divorce/input/
tail -n +2 divorce_phen.txt > divorce_phen.gcta

# get covs
cp -r /tsd/p805/data/durable/projects/Divorce_GWAS/Contin_cov_gcta.txt /cluster/projects/p805/rosac/divorce/input/
cp -r /tsd/p805/data/durable/projects/Divorce_GWAS/Catego_cov_gcta.txt /cluster/projects/p805/rosac/divorce/input/


# make ID keep file for grm  
cd /cluster/projects/p805/rosac/divorce/input/
awk '{ print $1, $2 }' divorce_phen.txt | tail -n +2 > divorce_phen.keep

# make a GRM to be used as a covariate in the GWAS, to control for relatedness.

/cluster/projects/p805/rosac/software/gcta_1.91.7beta/gcta64 \
--bfile /cluster/projects/p805/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc_regenie_500k_snps \
--make-grm \
--keep divorce_phen.keep \
--maf 0.01 \
--threads 10 \
--out Divorce_GRM


# make sparse GRM 
/cluster/projects/p805/rosac/software/gcta_1.91.7beta/gcta64 \
--grm Divorce_GRM \
--make-bK 0.05 \
--out Divorce_GRMK


# main GWA analysis in fastGWA controlling for GRM:
inp="/cluster/projects/p805/rosac/divorce/input/"
out="/cluster/projects/p805/rosac/divorce/output/"
cd ${inp}
/cluster/projects/p805/software/gcta/gcta-1.94.1 \
--fastGWA-mlm-binary \
--bfile /cluster/projects/p805/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
--grm-sparse Divorce_spGRM \
--pheno divorce_phen.gcta \
--mpheno 1 \
--qcovar Contin_cov_gcta.txt \
--covar Catego_cov_gcta.txt \
--est-vg HE \
--threads 10 \
--out ${out}Dissolve

# female
/cluster/projects/p805/software/gcta/gcta-1.94.1 \
--fastGWA-mlm-binary \
--bfile /cluster/projects/p805/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
--grm-sparse Divorce_spGRM \
--pheno female_divorce_pheno.gcta \
--mpheno 1 \
--qcovar Contin_cov_gcta.txt \
--covar Catego_cov_gcta.txt \
--est-vg HE \
--threads 10 \
--out ${out}Dissolve_Female

# male
/cluster/projects/p805/software/gcta/gcta-1.94.1 \
--fastGWA-mlm-binary \
--bfile /cluster/projects/p805/data/genetics/MoBaPsychGen_v1/MoBaPsychGen_v1-ec-eur-batch-basic-qc \
--grm-sparse Divorce_spGRM \
--pheno male_divorce_pheno.gcta \
--mpheno 1 \
--qcovar Contin_cov_gcta.txt \
--covar Catego_cov_gcta.txt \
--est-vg HE \
--threads 10 \
--out ${out}Dissolve_Male


# LDSC analyses using sumstats in R to get SNP h2 and rG between sexes.

library(GenomicSEM)
library(data.table)


munge(files=c('Dissolve.fastGWA', 'Dissolve_Male.fastGWA','Dissolve_Female.fastGWA'), hm3 = "w_hm3.noMHC.snplist",trait.names=c('Dissolution','Dissolution_Male','Dissolution_Female'),N=c(128476,52035,76441))

ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
# 
traits <- c("Dissolution_Male.sumstats.gz","Dissolution_Female.sumstats.gz")
sample.prev <- c(0.47,0.44)
population.prev <- c(0.392,0.392)


LDSCoutput <- ldsc(traits, sample.prev, population.prev, ld, wld)
