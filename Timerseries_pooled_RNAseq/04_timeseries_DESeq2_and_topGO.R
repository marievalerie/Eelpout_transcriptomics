library(DESeq2)
library(PCAtools)
library(apeglm)
library(tximport)
library(readr)
library(sva)
library(DEGreport)

#plotting 
library(devtools)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(grid)
library("pheatmap")
library("RColorBrewer")
library("ggcorrplot")



sig_lvl <- .05  #you can adjust this

#load data with tximport
#set to ballgown coverage data dir

setwd('../genome_guided_expression_gff_as_ref/')
files = c(read.table('samples_lst.txt')) ##file with .ctab file names
files_names <- c(files$V1)
files <- c(files$V2)
names(files) <- files_names

tx2gene <- read.table('/share/pool/mbrasseur/eelpout_RNA_timeseries/reference/zoarces_tx2gene')

#read in tx abundances; for timeseries -> 150 bp
txi.stringtie <- tximport(files, type = "stringtie", tx2gene = tx2gene, readLength = 150)

#change to analysis dir
setwd('/share/pool/mbrasseur/eelpout_RNA_timeseries/DESeq')

#load metadata
sampleData <- read.csv('timeseries_metadata.csv', sep = ';', header = 1)
row.names(sampleData) <- sampleData$Sample
sampleData <- sampleData[colnames(txi.stringtie$counts),]

sampleData$Year <- factor(sampleData$Year, levels = c(1994, 1999, 2005,2010,2012,2017,2019,2021)) ###treatment should be factor
levels(sampleData$Year) ##check the levels -> control must come first

sampleData$Location <- factor(sampleData$Location, levels = c("Dar", "Mel", "VaM"))

#check; if TRUE, proceed
all(colnames(txi.stringtie$counts) == rownames(sampleData)) 

#create DESeq object
dds1 <- DESeqDataSetFromTximport(txi.stringtie, colData=sampleData, design= ~ Location + Year) ##could make sense to fit time as numerical predictor an the add interaction term 

#to remove the mito rRNA/tRNA genes
#dds1 <- dds1[setdiff(rownames(dds1), '.'),]


#some filtering steps
dds1 <- estimateSizeFactors(dds1)
nc <- counts(dds1, normalized=TRUE)

filter2 <- rowSums(nc >= 3) >= 8 
dds1 <- dds1[filter2,]

#how many genes
length(row.names(counts(dds1))) ## 20041 genes 


####run DESeq() on the complete data set and use the LRT to test for the effect of time
dds1 <- DESeq(dds1, test="LRT", reduced = ~ Location) ##this will test if all time factor levels significantly contribute to the likelihood of the mdodel


#check dispersion
pdf(plot, file="gene_dispersion.pdf")
plot <- plotDispEsts(dds1)
dev.off()

#define color and shapes  for PCA
cols = c('#4B006E', 'darkgreen', 'darkorange')
shps = c(19,15)

#apply variance stabilzation on count data for PCA
vsd <- vst(dds1)

####PCA on adjusted counts
p <- pca(assay(vsd), metadata = colData(dds1))

pdf('PCA.pdf', height = 7, width = 5)
biplot(p, legendPosition = 'bottom', colby = 'Location', colkey = cols, shape = 'Sequence_batch', shapekey = shps, shapeLegendTitle = 'Sequencing batch',
       lab = p$metadata$Year, 
       xlab = paste0("PC1", ", ", round(p$variance["PC1"], digits = 2), "% variance"),
       ylab = paste0("PC2", ", ", round(p$variance["PC2"], digits = 2), "% variance")) 
dev.off()




####PCA reveals batch effect -> control for it
####surrogate variable analysis (SVA) to estimate latent factors

#extract normalized counts from DESeq object
dat  <- counts(dds1, normalized = TRUE)

#estimate latent factors i.e., surrogates
mod  <- model.matrix(~ Location + Year, colData(dds1)) ##specify full model
mod0 <- model.matrix(~   1, colData(dds1)) ##specifiy reduced model

#run sva
svseq <- svaseq(dat, mod, mod0) ##note down here the number of surrogates; 4


#frozen SVA to regress out latent factors 
newV = NULL ## neccessary due to bug in the sva pacakge
dat_vst  <- vst(dds1)
fsvaobj <- fsva(assay(dat_vst), mod, svseq, method = 'exact') 

#get the adjusted data for PCA to check
data_adj <- fsvaobj$db

#run pca again
p <- pca(data_adj, metadata = sampleData)

pdf('PCA_SVA_corrected.pdf', height = 7, width = 6)
biplot(p, legendPosition = 'bottom', colby = 'Location', colkey = cols, 
       lab = p$metadata$Year, 
       xlab = paste0("PC1", ", ", round(p$variance["PC1"], digits = 2), "% variance"),
       ylab = paste0("PC2", ", ", round(p$variance["PC2"], digits = 2), "% variance")) 
dev.off()


#### add SVs to DESeq2 model
ddssva <- dds1
ddssva$SV1 <- svseq$sv[,1] ##SV1
ddssva$SV2 <- svseq$sv[,2] ##SV2
ddssva$SV3 <- svseq$sv[,3] ##SV3
ddssva$SV4 <- svseq$sv[,4] ##SV4

design(ddssva) <- ~ SV1+SV2+SV3+SV4+Location+Year
ddssva <- DESeq(ddssva, test="LRT", reduced = ~ SV1+SV2+SV3+SV4+Location)

#remove the 30 genes that did not converge
ddsClean <- ddssva[which(mcols(ddssva)$betaConv),]
ddsClean <- ddssva[which(mcols(ddssva)$fullBetaConv),]

####batch effect corrected results
res_time_LRT <- results(ddsClean, alpha=sig_lvl)
summary(res_time_LRT)
'''
out of 20011 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 158, 0.79%
LFC < 0 (down)     : 358, 1.8%
outliers [1]       : 0, 0%
low counts [2]     : 2716, 14%
'''
write.csv(as.data.frame(res_time_LRT), file= "DEG_zoarces_timeseries_LRT.csv")




####test for DE between location with Wald test
ddssva2 <- DESeq(ddssva)
ddsClean2 <- ddssva2[which(mcols(ddssva2)$betaConv),]


#Varel-Mellum vs Darss
res.VaMvsDar <- results(ddsClean2, name = "Location_VaM_vs_Dar", alpha=sig_lvl)
'''
> summary(res.VaMvsDar)

out of 20011 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 589, 2.9%
LFC < 0 (down)     : 545, 2.7%
outliers [1]       : 0, 0%
low counts [2]     : 3104, 16%
'''

#apply effect size shrinkage
res.VaMvsDar.apeglm <- lfcShrink(dds=ddsClean2, type="apeglm", coef=7, res = res.VaMvsDar) 
write.csv(as.data.frame(res.VaMvsDar.apeglm), file= "DEG_zoarces_VaMvsDar_apeglm.csv")


#Meldorf Bay vs Darss
res.MelvsDar <- results(ddsClean2, name = "Location_Mel_vs_Dar", alpha=sig_lvl)
'''
out of 20011 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1568, 7.8%
LFC < 0 (down)     : 1384, 6.9%
outliers [1]       : 0, 0%
low counts [2]     : 1164, 5.8%
'''

#apply effect size shrinkage
res.MelvsDar.apeglm <- lfcShrink(dds=ddsClean2, type="apeglm", coef=6, res = res.MelvsDar) 
write.csv(as.data.frame(res.MelvsDar.apeglm), file= "DEG_zoarces_MelvsDar_apeglm.csv")







####clustering of time responsive genes according to their expression behaviour 
#extract the genes identified with the LRT in DESeq2
time_genes <- res_time_LRT
time_genes$padj[is.na(res_time_LRT$padj)] <- 1
time_genes <- time_genes[time_genes$padj < sig_lvl,]

#cluster based on batch effect adjusted counts
vsd_sig <- data_adj[row.names(time_genes),]

#preliminary
pdf('time_clusters.pdf')
cluster <- degPatterns(vsd_sig, sampleData, time = "Year", col = 'Marine_region', minc = 20)
dev.off()

#get the data to make a pretty plot
clustering_data <- cluster[["normalized"]]

#prepare labeller for plotting with facet wrap;
#you can count the number of observations with dplyr pipe but the observatiosn are not unique in the long format

clustering_data  %>% count(cluster) %>%
'''
cluster    n
1        1  576 
2        4 1280
3        5  352
4        7  480
5        8  384
6       12  400
7       13  688
8       14  416
9       15  416
10      16  432
'''

#each gene has 16 observations because odf 8 Northern Sea observations (16 Northern Sea samples were merged) and 8 Baltic Sea observations;
#therefore, these numbers must be divided by 16 to get the number of genes in each cluster

clustering_data  %>% count(cluster) %>%
  {. ->> geneCluster } 

geneCluster$n <- geneCluster$n/16

#define label function 

label_facet <- function(original_var, custom_name){
  lev <- levels(as.factor(original_var))
  lab <- paste0("Cluster ", original_var,": ", custom_name, ' genes')
  names(lab) <- lev
  return(lab)  
}

label_facet(geneCluster$cluster, geneCluster$n)


#connect mean z-values to get sinmething like a discrete trendline

pdf('time_cluster_with_discrete_trendline_both_regions.pdf')
ggplot(clustering_data,
       aes(Year, value, color = Marine_region)) +
  #geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  geom_boxplot(aes(Year, value)) +
  scale_color_manual(name = "Marine region", labels = c("Baltic Sea", "North Sea"), values = c('#4B006E', 'darkorange')) +
  labs(y= "z-score of expression", x = "") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_summary(aes(x=as.numeric(Year), color= Marine_region), fun=median, geom='line', linewidth=.8) +
  facet_wrap(vars(cluster), labeller = labeller(cluster = label_facet(geneCluster$cluster, geneCluster$n)), ncol = 2)
dev.off()


#find genes in cluster 14 (old name) / 8 (newname)
unique(clustering_data[clustering_data$cluster == 14, ]$genes)

'''
 [1] "g10860" "g1096"  "g11687" "g12955" "g13687" "g1375"  "g1645"  "g17031"
 [9] "g17333" "g22326" "g23107" "g23235" "g23346" "g2403"  "g2476"  "g26200"
[17] "g26399" "g27198" "g4985"  "g5031"  "g5095"  "g5111"  "g5225"  "g6375" 
[25] "g7048"  "g9660" 
'''


#find genes in cluster 7(old)/4(new)
unique(clustering_data[clustering_data$cluster == 7, ]$genes)

'''
 [1] "g10041" "g10609" "g1082"  "g12394" "g12798" "g15604" "g18501" "g18503"
 [9] "g23363" "g23555" "g24778" "g24828" "g2614"  "g27000" "g27450" "g3112" 
[17] "g3116"  "g3882"  "g4400"  "g5570"  "g5861"  "g6020"  "g7123"  "g7641" 
[25] "g8197"  "g8550"  "g8661"  "g9245"  "g9659"  "g9863" 
'''

#g10588, g11539, g14465, g16053, g17213, g6376 are solute carriers



####functional enrichment with topGO
library(tidyr)
library(GO.db)
library(topGO)

#first: test the timeseries genes against the DE background
#gene2go must be custom prepared prepared;transform to longformat

GO_ids = read.csv('StringtieEggNog_Gene2GO.csv', sep=';', header = F)
long_GO <- gather(GO_ids, Gen, IDs, V2:V1160)

#take out genes without GO terms
long_GO <- long_GO[which(long_GO$`IDs` != ""),] 

#remove variable column
long_GO <- long_GO[, c(1, 3)]

#sort by transcript/gene
gene.go <- long_GO[order(long_GO$V1), ]

#Create list with element for each gene, containing vectors with all terms for each gene
gene2GO <- tapply(gene.go$`IDs`, gene.go$V1, function(x)x)

head(gene2GO) #go IDs as strings

#extract time responsive genes
DE <- res_time_LRT
DE$padj[is.na(DE$padj)] <- 1

#define alpha for functional enrichment
pcutoff = 0.05 

DE$sig <- ifelse(DE$padj < pcutoff, 1, 0)
geneList <- DE$sig

#geneList need the same names as in the match for GO terms (gene2GO)
names(geneList) <- row.names(DE)

####run this code for the different ontologies you are interested in
                  
#Create topGOdata object:
GOdata <- new("topGOdata",
                   ontology = "BP",  #ontology criteria = biological process
                   allGenes = geneList, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                   geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

#run enrichment test (here Fishers Exact Test) and explore different algorithms to decorrelate the GO graph structure
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultFisher.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

#summarize in table
GOtable <- GenTable(GOdata, p.value.Fisher.elim = resultFisher.elim, p.value.Fisher.weight01 = resultFisher.weight01, topNodes = length(resultFisher.weight01@score),numChar = 10000000)

write.csv(GOtable, file = 'DESeq_timeseries_BP_enrichment_fisher.csv')
#write.csv(GOtable, file = 'DESeq_timeseries_MF_enrichment_fisher.csv')




####functional enrichment for DEG due to location 
#results are stored in objects res.VaMvsDar and res.MelvsDar; run the code for both

DE <- res.MelvsDar  
head(DE)
DE$padj[is.na(DE$padj)] <- 1

#differentiate between upregulated and downregulated to test the seperately
DE_up <- DE
DE_up$up <-ifelse(DE$log2FoldChange < 0, 0, 1)
DE_up$padj <- ifelse(DE_up$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_up$padj == 1 & DE_up$up == 1, 1, 0)
geneList_up <- tmp

DE_down <- DE
DE_down$down <-ifelse(DE_down$log2FoldChange > 0, 0, 1)
DE_down$padj <- ifelse(DE_down$padj < pcutoff, 1, 0)
tmp <- ifelse(DE_down$padj == 1 & DE_down$down == 1, 1, 0)
geneList_down <- tmp


#geneList need the same names as in the match for GO terms (gene2GO)
names(geneList_up) <- unlist(lapply(row.names(DE_up), function(x)x[1]))
names(geneList_down) <- unlist(lapply(row.names(DE_down), function(x)x[1]))

#Create topGOdata object:
GOdata_down <- new("topGOdata",
                   ontology = "MF",  #ontology criteria 
                   allGenes = geneList_down, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                   geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

#run enrichment test (here Fishers Exact Test) and explore different algorithms to decorrelate the GO graph structure
resultFisher_down.elim <- runTest(GOdata_down, algorithm = "elim", statistic = "fisher")
resultFisher_down.weight01 <- runTest(GOdata_down, algorithm = "weight01", statistic = "fisher")


#summarize in table
down <- GenTable(GOdata_down, p.value.elim = resultFisher_down.elim, p.value.weight01 = resultFisher_down.weight01, topNodes = length(resultFisher_down.elim@score), numChar = 1000000)
down$regulation <- 'down'


GOdata_up <- new("topGOdata",
                 ontology = "MF",  #ontology criteria
                 allGenes = geneList_up, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                 geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                 annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms


resultFisher_up.elim <- runTest(GOdata_up, algorithm = "elim", statistic = "fisher")
resultFisher_up.weight01 <- runTest(GOdata_up, algorithm = "weight01", statistic = "fisher")


up <- GenTable(GOdata_up, p.value.elim = resultFisher_up.elim, p.value.weight01 = resultFisher_up.weight01, topNodes = length(resultFisher_up.elim@score), numChar = 10000000)
up$regulation <- 'up'

#combine results
results <- rbind (down, up)
results <- results[results$Annotated > 3,]

write.csv(results, file = 'MelvsDar_MF_enrichment.csv')


#find the genes behind 3 putative salinity go terms

GOdata_up <- new("topGOdata",
                 ontology = "BP",  #ontology criteria
                 allGenes = geneList_up, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                 geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                 annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

selected <- c("GO:0098656", "GO:1903276", "GO:0033561") 
sigGenes(GOdata_up)
sigGenes(GOdata_up)[sigGenes(GOdata_up)  %in% unlist(unname(genesInTerm(GOdata_up, selected[1])))] ##the order matters here
#g10588, g11539, g14465, g16053, g17213, g6376 are solute carriers


##############################################

####check, how the heat responsive genes behave in the timeseries data
#load the data
heatstress <- read.csv('../../eelpout_RNA_heatstress_experiment/DESeq/DEG_zoarces_heatstress_exp.csv', row.names = 1)

heatstress$padj[is.na(heatstress$padj)] <- 1
heatstress <- heatstress[heatstress$padj < 0.05,]
heatstress <- row.names(heatstress)

heatstress <- intersect(heatstress, row.names(data_adj)) ##2440 genes
vsd_heatstress <- data_adj[heatstress,]

#preliminary clustering
pdf('heatstress_clusters_in_timeseries_data.pdf')
cluster <- degPatterns(vsd_heatstress, sampleData, time = "Year", col = 'Marine_region', minc = 20)
dev.off()

#get the data                                      
clustering_data <- cluster[["normalized"]]
clustering_data  %>% count(cluster) 

'''
   cluster    n
1        2  464
2        3 1744
3        4 1120
4        5 2176
5        6  464
6        9 1392
7       10  640
8       11 1792
9       12  384
10      15 1216
11      19  368
12      20  688
13      22  672
14      23  512
15      24  384
16      25  784
17      26 1088
18      29  576
19      31  336
20      32  784
21      37  368
22      38  336
23      40  592
24      42  976
25      45  336
26      51  368
27      53  640
28      58  464
29      64  416
30      65  416
31      80  400
32      95  336
33     103  336
34     112  432
'''


clustering_data  %>% count(cluster) %>%
  {. ->> geneCluster } 

geneCluster$n <- geneCluster$n/16

#order by size to show only largest cluster profiles in the MS
geneCluster <- geneCluster[order(geneCluster$n, decreasing=T),]

#large clusters (here:more than 40 genes)
geneCluster.large <- geneCluster[geneCluster$n >= 40,]
clustering_data.large <- clustering_data[clustering_data$cluster %in% geneCluster.large$cluster,]

#small clusters
geneCluster.small <- geneCluster[geneCluster$n < 40,]
clustering_data.small <- clustering_data[clustering_data$cluster %in% geneCluster.small$cluster,]


#define labels 
label_facet(geneCluster$cluster, geneCluster$n)

#large clusters
pdf('heatstress_cluster_in_timeseries_with_discrete_trendline_both_regions_large.pdf')
ggplot(clustering_data.large,
       aes(Year, value, color = Marine_region)) +
  #geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  geom_boxplot(aes(Year, value)) +
  scale_color_manual(name = "Marine region", labels = c("Baltic Sea", "North Sea"), values = c('#4B006E', 'darkorange')) +
  labs(y= "z-score of expression", x = "") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_summary(aes(x=as.numeric(Year), color= Marine_region), fun=median, geom='line', linewidth=.8) +
  facet_wrap(vars(cluster), ncol = 4,
             labeller = labeller(cluster = label_facet(geneCluster.large$cluster, geneCluster.large$n)))
dev.off()


geneCluster.small <- geneCluster[geneCluster$n < 40,]
clustering_data.small <- clustering_data[clustering_data$cluster %in% geneCluster.small$cluster,]

#small clusters
pdf('heatstress_cluster_in_timeseries_with_discrete_trendline_both_regions_small.pdf')
ggplot(clustering_data.small,
       aes(Year, value, color = Marine_region)) +
  #geom_point(position = position_jitterdodge(dodge.width = 0.9)) +
  geom_boxplot(aes(Year, value)) +
  scale_color_manual(name = "Marine region", labels = c("Baltic Sea", "North Sea"), values = c('#4B006E', 'darkorange')) +
  labs(y= "z-score of expression", x = "") +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_summary(aes(x=as.numeric(Year), color= Marine_region), fun=median, geom='line', linewidth=.8) +
  facet_wrap(vars(cluster), ncol = 4,
             labeller = labeller(cluster = label_facet(geneCluster.small$cluster, geneCluster.small$n)))
dev.off()
