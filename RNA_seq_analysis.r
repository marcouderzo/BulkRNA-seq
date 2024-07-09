# R Code

library("Rsubread")
setwd("/home/ubuntu/data/mydatalocal/TP_RNASeq_genome/")
bam_files <- c("FL_ko_1.bam", "FL_ko_2.bam", "FL_wt_1.bam", "FL_wt_2.bam",
               "T3_ko_1.bam", "T3_ko_2.bam", "T3_wt_1.bam", "T3_wt_2.bam")
bam_files <- system("ls align_out/*bam", intern=TRUE)
gtf_file <- "annot/annot.gtf"
counts <- featureCounts(files = bam_files, annot.ext = gtf_file, isGTFAnnotationFile = TRUE)
counts$counts["ENSMUSG00000086903",]

#Out:
#FL_ko_1.bam FL_ko_2.bam FL_wt_1.bam FL_wt_2.bam T3_ko_1.bam T3_ko_2.bam
#          0           0           5           8           0           0
#T3_wt_1.bam T3_wt_2.bam
#        427         313

#Can confirm the deletion of the Hotair gene in knockout compared 
#to wildtype, and it is much more expressed in trunk rather than in four limbs.

length <- counts$annotation$Length
ex_length_kb <- counts$annotation$Length / 10^-3 # I have one exonic length per gene
sum_RPKs <- colSums(counts$counts)
cnts <- counts$counts
RPK <- cnts/ex_length_kb
TPM <- t(RPK)/(sum_RPKs)/(10^6)


# Differential Expression Analysis

read.counts <- read.table("ReadCounts_AllGenes.txt")
print(head(read.counts))
                                     GeneID GeneName ExonicLength FL_wt_1
#ENSMUSG00000000001_Gnai3 ENSMUSG00000000001    Gnai3         3262    8715
#ENSMUSG00000000003_Pbsn  ENSMUSG00000000003     Pbsn          902       0
#ENSMUSG00000000028_Cdc45 ENSMUSG00000000028    Cdc45         2252    1862
#ENSMUSG00000000031_H19   ENSMUSG00000000031      H19         2372  210565
#ENSMUSG00000000037_Scml2 ENSMUSG00000000037    Scml2         6039     224
#ENSMUSG00000000049_Apoh  ENSMUSG00000000049     Apoh         1594       0
#                         FL_wt_2 FL_ko_1 FL_ko_2 T3_wt_1 T3_wt_2 T3_ko_1
#ENSMUSG00000000001_Gnai3    9466    7800   10837    9605    8742   10144
#ENSMUSG00000000003_Pbsn        0       0       0       0       0       0
#ENSMUSG00000000028_Cdc45    1826    1929    2453    1717    1481    1738
#ENSMUSG00000000031_H19    178310  152329  223100  129665  125728   85720
#ENSMUSG00000000037_Scml2     263     195     265     237     196     164
#ENSMUSG00000000049_Apoh        1       0       0       0       0       0
#                         T3_ko_2
#ENSMUSG00000000001_Gnai3   11636
#ENSMUSG00000000003_Pbsn        0
#ENSMUSG00000000028_Cdc45    2120
#ENSMUSG00000000031_H19     99315
#ENSMUSG00000000037_Scml2     160
#ENSMUSG00000000049_Apoh        0


read.counts.t3 <- read.counts[, grepl("^T3_", names(read.counts))]
print(head(read.counts.t3))


#                         T3_wt_1 T3_wt_2 T3_ko_1 T3_ko_2
#ENSMUSG00000000001_Gnai3    9605    8742   10144   11636
#ENSMUSG00000000003_Pbsn        0       0       0       0
#ENSMUSG00000000028_Cdc45    1717    1481    1738    2120
#ENSMUSG00000000031_H19    129665  125728   85720   99315
#ENSMUSG00000000037_Scml2     237     196     164     160
#ENSMUSG00000000049_Apoh        0       0       0       0
coldata.t3 <- data.frame("genotype")
View(coldata.t3)
genotype <- c("wt", "ko")
coldata.t3 <- data.frame()
genotype <- c("wt", "wt", "ko", "ko")
coldata.t3 = data.frame(genotype=genotype)


print(coldata.t3)
#  genotype
#1       wt
#2       wt
#3       ko
#4       ko
#> individual.t3 = samplw_names[8:11]
#Erreur : objet 'samplw_names' introuvable
#> individual.t3 = sample_names[8:11]
#Erreur : objet 'sample_names' introuvable
sample_names = colnames(read.counts)
rownames(coldata.t3)=sample_names[8:11]
coldata.t3

#        genotype
#T3_wt_1       wt
#T3_wt_2       wt
#T3_ko_1       ko
#T3_ko_2       ko

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = read.counts.t3, colData = coldata.t3, design= ~ genotype)
dds <- DESeq(dds)
res <- results(dds)
res

#estimating size factors
#estimating dispersions
#gene-wise dispersion estimates
#mean-dispersion relationship
#final dispersion estimates
#fitting model and testing
#log2 fold change (MLE): genotype wt vs ko
#Wald test p-value: genotype wt vs ko
#DataFrame with 46425 rows and 6 columns
#                                   baseMean log2FoldChange     lfcSE      stat
#                                  <numeric>      <numeric> <numeric> <numeric>
#ENSMUSG00000000001_Gnai3          10012.693      -0.267351 0.1037334  -2.57729
#ENSMUSG00000000003_Pbsn               0.000             NA        NA        NA
#ENSMUSG00000000028_Cdc45           1759.199      -0.289040 0.1342800  -2.15251
#ENSMUSG00000000031_H19           109719.708       0.445782 0.0959334   4.64678
#ENSMUSG00000000037_Scml2            188.848       0.393353 0.2694431   1.45987
#...                                     ...            ...       ...       ...
#ENSMUSG00000108295_RP24-434O21.2       0.00             NA        NA        NA
#ENSMUSG00000108296_RP23-373P9.1        0.00             NA        NA        NA
#ENSMUSG00000108297_RP23-392K24.6       0.75        3.01728   4.47387  0.674424
#ENSMUSG00000108298_RP24-315D19.3       0.00             NA        NA        NA
#ENSMUSG00000108299_RP23-120I17.1       0.00             NA        NA        NA
#                                      pvalue        padj
#                                   <numeric>   <numeric>
#ENSMUSG00000000001_Gnai3         9.95779e-03 0.055597977
#ENSMUSG00000000003_Pbsn                   NA          NA
#ENSMUSG00000000028_Cdc45         3.13569e-02 0.125227628
#ENSMUSG00000000031_H19           3.37150e-06 0.000123043
#ENSMUSG00000000037_Scml2         1.44325e-01 0.341205260
#...                                      ...         ...
#ENSMUSG00000108295_RP24-434O21.2          NA          NA
#ENSMUSG00000108296_RP23-373P9.1           NA          NA
#ENSMUSG00000108297_RP23-392K24.6    0.500042          NA
#ENSMUSG00000108298_RP24-315D19.3          NA          NA
#

res_filtered <- res[!is.na(res$padj), ]
alpha <- 0.05
significant_genes <- sum(res_filtered$padj < alpha)
significant_genes

#2574 genes, based on adjusted p-value < alpha.

plot <- plotMA(res, ylim=c(-2,2))
resultsNames(dds)
#[1] "Intercept"         "genotype_wt_vs_ko"
resLFC <- lfcShrink(dds, coef="genotype_wt_vs_ko", type="apeglm")

# Same thing for FL

read.counts.fl <- read.counts[, grepl("^FL_", names(read.counts))]
print(head(read.counts.fl))


coldata.fl <- data.frame("genotype")
View(coldata.fl)
genotype <- c("wt", "ko")
coldata.fl <- data.frame()
genotype <- c("wt", "wt", "ko", "ko")
coldata.fl = data.frame(genotype=genotype)


print(coldata.fl)

sample_names = colnames(read.counts)
rownames(coldata.fl)=sample_names[4:7]
coldata <- t(coldata.fl)


library("DESeq2")
dds_fl <- DESeqDataSetFromMatrix(countData = read.counts.fl, colData = coldata.fl, design= ~ genotype)
dds_fl <- DESeq(dds_fl)
res_fl <- results(dds_fl)
res_fl

res_fl_filtered <- res_fl[!is.na(res_fl$padj), ]
alpha <- 0.05
significant_genes <- sum(res_fl_filtered$padj < alpha)
significant_genes


plot_fl <- plotMA(res_fl, ylim=c(-2,2))


resultsNames(dds_fl)
#[1] "Intercept"         "genotype_wt_vs_ko"
res_flLFC <- lfcShrink(dds_fl, coef="genotype_wt_vs_ko", type="apeglm")