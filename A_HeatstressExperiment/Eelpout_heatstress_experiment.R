#analysis
library(DESeq2)
library(PCAtools)
library(apeglm)
library(tximport)
library(readr)
library(edgeR)
library(limma)
#plotting 
library(devtools)
library(reshape2)
library(ggplot2)
#library(tidyverse)
library(grid)
library(dplyr)
library("pheatmap")
library("RColorBrewer")
library("ggcorrplot")

sig_lvl <- .05  ##you can adjust this

###load data with tximport

#set to ballgown coverage data dir
setwd('../genome_guided_expression_gff_as_ref/')
files = c(read.table('samples.txt')) ##file with .ctab file names
files_names <- c(files$V1)
files <- c(files$V2)
names(files) <- files_names

#create tx2gene file
a <- read_tsv(files[1]) ##somehow does not work when I want to store it in temp :/
tx2gene <- a[, c("t_name", "gene_id")]
tx2gene <- tx2gene[tx2gene$gene_id != ".",] ##check the mitogene issue!!!

#thisdoes not contain mitogenes
#tx2gene <- read.table('/share/pool/mbrasseur/eelpout_RNA_timeseries/reference/zoarces_tx2gene')

##read in tx abundances for deseq
txi.stringtie <- tximport(files, type = "stringtie", tx2gene = tx2gene, readLength = 100)

#change to analysis dir
setwd('/share/pool/mbrasseur/eelpout_RNA_heatstress_experiment/DESeq')

##load metadata
sampleData <- read.csv('eelpout_heatstress_exp_treatments.csv', sep = ';', header = 1)
row.names(sampleData) <- sampleData$Library_name
sampleData <- sampleData[colnames(txi.stringtie$counts),]

sampleData$Temperature_treatment <- factor(sampleData$Temperature_treatment) ###treatment should be factor
levels(sampleData$Temperature_treatment) ##check the levels -> control must come first


##check; if TRUE, proceed
all(colnames(txi.stringtie$counts) == rownames(sampleData)) 

#runDESeq
dds1 <- DESeqDataSetFromTximport(txi.stringtie, colData=sampleData, design= ~ Temperature_treatment)

##remove the mito rRNA/tRNA genes
#dds1 <- dds1[setdiff(rownames(dds1), '.'),]


##some filtering steps; omit this for now
dds1 <- estimateSizeFactors(dds1)
nc <- counts(dds1, normalized=TRUE)

filter2 <- rowSums(nc >= 3) >= 6 
dds1 <- dds1[filter2,]

##no. of genes (after filtering)
length(row.names(counts(dds1))) ## 17411 genes


##run DESeq() on the complete data set
dds1 <- DESeq(dds1) 

##check dispersion
pdf(plot, file="gene_dispersion.pdf")
plot <- plotDispEsts(dds1)
dev.off()


### define some colors for the PCA to differentiate between colors
cols = c('skyblue3',  'palevioletred3')

##apply variance stabilzation on count data for PCA
vsd <- vst(dds1)

##pca
p <- pca(assay(vsd), metadata = colData(dds1))

##plot the first 2 axes in a biplot (if you want to hide sample labels, add lab = NULL)
pdf('PCA.pdf', height = 7, width = 5)
biplot(p, legendPosition = 'bottom', colby = 'Temperature_treatment', colkey = c('black', 'black'), shape = 'Temperature_treatment', shapekey = c(0,7), shapeLegendTitle = 'Temperature',
       lab = p$metadata$Tank, 
       xlab = paste0("PC1", ", ", round(p$variance["PC1"], digits = 2), "% variance"),
       ylab = paste0("PC2", ", ", round(p$variance["PC2"], digits = 2), "% variance")) 
dev.off()


##sample clustering based on distances
sampleDists <- dist(t(assay(vsd)), method = 'man')
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Library_name, vsd$Tank, vsd$Temperature_treatment, sep = "_" )
colnames(sampleDistMatrix) <- NULL
svg('manhatten_distance_heatmap.svg', height = 6, width = 7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255))
dev.off()

###get the results
res <- results(dds1, name = "Temperature_treatment_increased_vs_control", lfcThreshold = 0, alpha = sig_lvl)

##you can get a summary of the results:
summary(res)

'''
out of 17412 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1198, 6.9%
LFC < 0 (down)     : 1243, 7.1%
outliers [1]       : 86, 0.49%
low counts [2]     : 1350, 7.8%
'''


##save this
write.csv(as.data.frame(res), file= "DEG_zoarces_heatstress_exp.csv") 

##shrink the LFC (wich are effect sizes)
res_apeglm <- lfcShrink(dds=dds1, type="apeglm", coef=2, res = res) 

#save this
write.csv(as.data.frame(res_apeglm), file= "DEG_zoarces_heatstress_exp_apeglm.csv") 


##-> use this geneset!

####### plotting heatmap of normalized counts of DE genes
library(reshape2)
library(viridis)
library(ComplexHeatmap)
library(dendsort)
library(circlize)

deseq2VST <- as.data.frame(assay(vsd))
deseq2VST$Gene <- rownames(deseq2VST)

##get significant genes

genes <- res
genes$padj[is.na(genes$padj)] <- 1
genes <- genes[genes$padj < sig_lvl,]

sig_genes <- row.names(genes)

deseq2VST <- deseq2VST[deseq2VST$Gene %in% sig_genes,]


# Make a heatmap
row_dend = dendsort(hclust(dist(as.matrix(deseq2VST[,c(1:12)]), method = 'man'))) 
col_dend = dendsort(hclust(dist(t(as.matrix(deseq2VST[,c(1:12)])), method = 'man')))

z.score <- t(scale(t(deseq2VST[,c(1:12)])))

col_labels <- c(paste('rep', c(1:6), sep = '. '), paste('rep', c(1:6), sep = '. '))
f1 = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))


heatmap <- Heatmap(as.matrix(z.score),col=f1, use_raster = F, cluster_rows = row_dend, cluster_columns = col_dend,
                   column_labels = col_labels, column_names_gp = gpar(fontsize = 10), column_names_rot = 45, 
                   column_title = c('12? C', '18? C'), column_title_side = "bottom", 
                   show_heatmap_legend = TRUE, show_row_names = FALSE, column_split = 2, 
                   heatmap_legend_param = (
                     list(title = "Expression", title_gp = gpar(fontsize = 13), legend_direction = "horizontal", col_fun = f1)))
heatmap

svg("heatmap_zscores_new_label.svg", height = 8, width = 10)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()        


#############Volcano plot#############

##read the data
conc <- as.data.frame(res_apeglm)
conc$padj[is.na(conc$padj)] <- 1

##define chr vector for color mapping
keyvals <- ifelse(
  conc$log2FoldChange < 0, 'skyblue3','palevioletred3')
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'palevioletred3'] <- 'upregulated'
names(keyvals)[keyvals == 'black'] <- 'n.s.'
names(keyvals)[keyvals == 'skyblue3'] <- 'downregulated'

##add expression info to df
conc <- conc %>% 
  mutate(
    Expression = case_when(log2FoldChange > 0 & padj < sig_lvl ~ "Upregulated\n(1,198)",
                           log2FoldChange < 0 & padj < sig_lvl ~ "Downregulated\n(1,243)",
                           TRUE ~ "n.s.\n(14,970)"))

conc$Expression <- factor(conc$Expression, levels = c("Downregulated\n(1,243)","Upregulated\n(1,198)",'n.s.\n(14,970)'))

##check if everything is right
conc %>% 
  count(Expression)

pdf('volcano_plot.pdf', width=6, height = 7)
p1 <- ggplot(conc, aes(log2FoldChange, -log(padj,10))) + # -log10 conversion  
  geom_point(aes(color = Expression), size = 1) + xlim(-11.5, 11.5)  + ylim (0, 27) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"padj"))+
  scale_color_manual(values = c("skyblue3","palevioletred3", "gray50" )) + theme_bw() +
  theme(legend.title = element_blank(), legend.text=element_text(size=16), legend.position = "bottom", axis.text=element_text(size=14), axis.title=element_text(size=16),
        panel.border = element_rect(linewidth = 1.5), panel.grid.minor = element_line(linewidth = 0.75), panel.grid.major = element_line(linewidth = 0.75))+
  guides(colour = guide_legend(override.aes = list(size=4)))
p1
dev.off()



####check with limma the tank issue and model tank as random effect

### with gene data set filtered and normalized as with DESeq2 before
x <- counts(dds1, normalized=TRUE)
dim(x) ## 17411 genes

v <- voom(x, design)
cor <- duplicateCorrelation(v, design, block = tank)
cor$consensus #0.1745016
#final model, accounting for tank as random effect
y <- voom(x, design, block = tank, correlation=cor$consensus)
##linear model fit
fit <- lmFit(y, design, block=tank, correlation=cor$consensus)
#computes moderated t-stats, F-stats 
fit <- eBayes(fit)

summary(decideTests(fit, p.value = 0.1))

'''
> summary(decideTests(fit, p.value = 0.1))
       (Intercept) treatmentincreased
Down          4069                255
NotSig        1855              15492
Up           11487               1664
'''

##get the genes

sig_genes_limma_tank <- topTable(fit, number=1919,sort.by="p")


# do not account for tank as random effect.
y      <- voom(x,design)
fit    <- lmFit(y,design)
fit    <- eBayes(fit)

summary(decideTests(fit, p.value = 0.1))

'''
Down          4130                404
NotSig        1721              14566
Up           11560               2441
'''

sig_genes_limma <- topTable(fit, number=2845,sort.by="p")

library(VennDiagram)

sig_genes #DESe2
row.names(sig_genes_limma) #limma without tank as random effect
row.names(sig_genes_limma_tank) #limma with random effect

sets <- list(sig_genes, row.names(sig_genes_limma), row.names(sig_genes_limma_tank))

# Chart
venn.diagram(
  x = sets,
  category.names = c("DESeq2" , "limma: - tank effect" , "limma: + tank effect"),
  filename='venn_diagram.png',
  imagetype = 'png',
  #circles
  fill = c('#5D2586', '#F98128', '#50B848'),    
  lwd = 1,
  lty = 'blank',
  euler.d = T,
  scaled = T,
  #numbers
  cex = .7,
  fontface = "bold",
  fontfamily = "sans",
  print.mode = 'raw',
  #labels
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.075),
  rotation = 1)




##############################
###with edgeR normalization
#read in tx abundances
##for preliminary edgeR/limma voom integration
'''
txi.stringtie <- tximport(files, type = "stringtie", tx2gene = tx2gene, readLength = 100, countsFromAbundance = "lengthScaledTPM")
gene.counts <- txi.stringtie$counts

##edgeR normalization
dge <- DGEList(counts=gene.counts, group = sampleData$Temperature_treatment)
x <- dge
dim(x) #28023 genes
keep.exprs <- filterByExpr(x, group=sampleData$Temperature_treatment, min.count = 3)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) #17357 genes

#calculate size factors which can be used by voom (?)
x <- calcNormFactors(x, method = "TMM")

##get design for limma
treatment = sampleData[, 4]
tank = sampleData[, 2]

design <- model.matrix(~treatment)

##
v <- voom(x, design)
cor <- duplicateCorrelation(v, design, block = tank)
cor$consensus #0.1884557
#final model, accounting for tank as random effect
y <- voom(x, design, block = tank, correlation=cor$consensus)
##linear model fit
fit <- lmFit(y, design, block=tank, correlation=cor$consensus)
#computes moderated t-stats, F-stats 
fit <- eBayes(fit)
summary(decideTests(fit))
'''
'''
> summary(decideTests(fit))
       (Intercept) treatmentincreased
Down          3460                169
NotSig        2370              17098
Up           11527                 90
'''
'''
topTable(fit, number=260,sort.by="p")

# do not account for tank as random effect.
y      <- voom(x,design)
fit    <- lmFit(y,)
fit    <- eBayes(fit)

summary(decideTests(fit))
'''
'''
> summary(decideTests(fit))
       (Intercept) treatmentincreased
Down          3592                390
NotSig        2158              16741
Up           11607                226
'''


#########################################







##################################################subset the data
##dds2 will be a four replicates data set, containing
#959-7, 959-11 (temp_treat) & 959-5, 959-1 (control)

#-> find indices for 959-8, 959-12, 959-6,959-2
colData(dds1)

dds2 <- dds1[,-c(3,5,9,11)]
#check
colData(dds1)

##run DESeq() again to re-estimate paramteres
dds2 <- DESeq(dds2)

vsd <- vst(dds2)

##pca
p <- pca(assay(vsd), metadata = colData(dds2))

##plot the first 2 axes in a biplot (if you want to hide sample labels, add lab = NULL)
pdf('PCA.pdf')
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'right', colby = 'Temperature_treatment', colkey = cols, shape = 'Tank', shapekey = shps) 
dev.off()


###get the results
res2 <- results(dds2, name = "Temperature_treatment_increased_vs_control", lfcThreshold = 0, alpha = sig_lvl)

##you can get a summary of the results:
summary(res2)

'''
out of 16282 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 745, 4.6%
LFC < 0 (down)     : 815, 5%
outliers [1]       : 105, 0.64%
low counts [2]     : 1260, 7.7%
'''


##save this
write.csv(as.data.frame(res2), file= "DEG_zoarces_heatstress_exp_subset1.csv") 

##shrink the LFC (wich are effect sizes)
res_apeglm2 <- lfcShrink(dds=dds2, type="apeglm", coef=2, res = res2) 

#save this
write.csv(as.data.frame(res_apeglm2), file= "DEG_zoarces_heatstress_exp_apeglm_subset1.csv") 



##subset 2
##dds3 will be a four replicates data set, containing
#959-8, 959-12, 959-6,959-2


#-> find indices for 959-7, 959-11, 959-5, 959-1 
colData(dds1)

dds3 <- dds1[,-c(10,2,8,4)]
#check
colData(dds3)

##run DESeq() again to re-estimate paramteres
dds3 <- DESeq(dds3)

vsd <- vst(dds3)

##pca
p <- pca(assay(vsd), metadata = colData(dds3))

##plot the first 2 axes in a biplot (if you want to hide sample labels, add lab = NULL)
pdf('PCA.pdf')
biplot(p, showLoadings = F, sizeLoadingsNames = 2, boxedLoadingsNames = F, legendPosition = 'right', colby = 'Temperature_treatment', colkey = cols, shape = 'Tank', shapekey = shps) 
dev.off()


###get the results
res3 <- results(dds3, name = "Temperature_treatment_increased_vs_control", lfcThreshold = 0, alpha = sig_lvl)

##you can get a summary of the results:
summary(res3)

'''
out of 16282 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 529, 3.2%
LFC < 0 (down)     : 630, 3.9%
outliers [1]       : 102, 0.63%
low counts [2]     : 1260, 7.7%
'''


##save this
write.csv(as.data.frame(res3), file= "DEG_zoarces_heatstress_exp_subset2.csv") 

##shrink the LFC (wich are effect sizes)
res_apeglm3 <- lfcShrink(dds=dds3, type="apeglm", coef=2, res = res3) 

#save this
write.csv(as.data.frame(res_apeglm3), file= "DEG_zoarces_heatstress_exp_apeglm_subset2.csv") 



#####################correlate dispersion values of significant genes
data1 <- as.data.frame(cbind(rownames(dds1), dispersions(dds1), dispersions(dds2)))
colnames(data1) <- c('genes', 'full_dataset', 'data_subset1')
row.names(data1) <- data1$genes

##first: give all genes the color NA 
data1$color <- NA


##then assign color to genes which are DE in the subset1
res2$padj[is.na(res2$padj)] <- 1
genesDE_sub <- res2[res2$padj < sig_lvl, ]

data1[row.names(genesDE_sub), 'color'] <- 'sub'

##subset to genes being DE in the full data set
res$padj[is.na(res$padj)] <- 1
genesDE_full<- res[res$padj < sig_lvl, ]

data1[row.names(genesDE_full), 'color'] <- 'full'

###then assign color to genes which are shared DE
shared <- intersect(row.names(genesDE_sub), row.names(genesDE_full))
data1[shared, 'color'] <- 'shared'


##omit NA genes
data1 <- na.omit(data1)
length(data1$genes) #2738

pdf('corrplot1.pdf')
ggplot(data1, aes(x=full_dataset, y=data_subset1, color=color)) + geom_point() +
  theme_minimal(base_size = 15)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())
dev.off()



##data subset2
data2 <- as.data.frame(cbind(rownames(dds1), dispersions(dds1), dispersions(dds3)))
colnames(data2) <- c('genes', 'full_dataset', 'data_subset2')
row.names(data2) <- data2$genes

##first: give all genes the color NA 
data2$color <- NA


##then assign color to genes which are DE in the subset2
res3$padj[is.na(res3$padj)] <- 1
genesDE_sub <- res3[res3$padj < sig_lvl, ]

data2[row.names(genesDE_sub), 'color'] <- 'sub'

##subset to genes being DE in the full data set
#res$padj[is.na(res$padj)] <- 1
#genesDE_full<- res[res$padj < sig_lvl, ]

data2[row.names(genesDE_full), 'color'] <- 'full'

###then assign color to genes which are shared DE
shared <- intersect(row.names(genesDE_sub), row.names(genesDE_full))
data2[shared, 'color'] <- 'shared'


##omit NA genes
data2 <- na.omit(data2)
length(data2$genes) #2553

pdf('corrplot2.pdf')
ggplot(data2, aes(x=full_dataset, y=data_subset2, color=color)) + geom_point() +
  theme_minimal(base_size = 15)+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank())
dev.off()

#################################################################################











#################functional enrichment:
library(tidyr)
library(GO.db)
library(topGO)

setwd('C:/Users/mbras/Desktop/eelpout_heatstress_experiment/functional annotation')


#######if gene2go is already prepared
##check that no invisible line breaks are in the file ;open the file one time a just save it
#transform to longformat

#to please the reviewer
GO_ids = read.csv('StringtieEggNog_Gene2GO.csv', sep=';', header = F)

long_GO <- gather(GO_ids, Gen, IDs, V2:V1160)


# take out genes without GO terms
long_GO <- long_GO[which(long_GO$`IDs` != ""),] 

##remove variable column
long_GO <- long_GO[, c(1, 3)]

#sort by transcript/gene
gene.go <- long_GO[order(long_GO$V1), ]

# Create list with element for each gene, containing vectors with all terms for each gene
gene2GO <- tapply(gene.go$`IDs`, gene.go$V1, function(x)x)

head(gene2GO) #go IDs as strings

##differentiate between up and downregulated genes
DE <- read.csv('../DESeq_for_MS/analysis with all samples/DEG_zoarces_heatstress_exp.csv')

#DE$X <- row.names(DE)
DE <- DE[, c(1,3,7)]
DE$padj[is.na(DE$padj)] <- 1

head(DE)

pcutoff = 0.05 ##you can adjust this

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


##geneList need the same names as in the match for GO terms (gene2GO)
names(geneList_up) <- unlist(lapply(DE_up$X, function(x)x[1]))
names(geneList_down) <- unlist(lapply(DE_down$X, function(x)x[1]))

##Create topGOdata object:

GOdata_down <- new("topGOdata",
                   ontology = "MF",  #ontology criteria 
                   allGenes = geneList_down, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                   geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

##run enrichment test: here Fishers Exact Test
resultFisher_down.elim <- runTest(GOdata_down, algorithm = "elim", statistic = "fisher")
resultFisher_down.weight01 <- runTest(GOdata_down, algorithm = "weight01", statistic = "fisher")


#summarize in table
down <- GenTable(GOdata_down, p.value.elim = resultFisher_down.elim, p.value.weight01 = resultFisher_down.weight01, topNodes = length(resultFisher_down.elim@score), numChar = 1000000)

down <- down[as.numeric(down$p.value.elim) < .05 | as.numeric(down$p.value.weight01) < .05 ,]

down$regulation <- 'down'
 

GOdata_up <- new("topGOdata",
                 ontology = "MF",  #ontology criteria
                 allGenes = geneList_up, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                 geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                 annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

resultFisher_up.elim <- runTest(GOdata_up, algorithm = "elim", statistic = "fisher")
resultFisher_up.weight01 <- runTest(GOdata_up, algorithm = "weight01", statistic = "fisher")


up <- GenTable(GOdata_up, p.value.elim = resultFisher_up.elim, p.value.weight01 = resultFisher_up.weight01, topNodes = length(resultFisher_up.elim@score), numChar = 10000000)

up <- up[as.numeric(up$p.value.elim) < .05 | as.numeric(up$p.value.weight01) < .05,]
up$regulation <- 'up'

results <- rbind (down, up)

results <- results[results$Annotated > 3,]

write.csv(results, file = 'MS_results/Zoarces_MF_enrichment.csv')




#############visualization for henrik
library("ggplot2")
library("tidyverse")
library("ggpubr")

ggdata <- read.csv ("/Users/mariebrasseur/Desktop/trier/eelpout_heatstress_experiment/functional annotation/MS_results/Zoarces_BP_enrichment.csv", sep =';', row.names = 1)

#upregulated genes
ggdata_up <- ggdata[ggdata$regulation == 'up',]

ggdata_up <- ggdata_up %>% 
  arrange(p.value.elim) %>%  # arrange in descending order
  slice(1:10) ##top10    

ggdata_up

ggdata_up$Term <- factor(ggdata_up$Term, levels = rev(ggdata_up$Term)) # fixes order
ggdata_up$p.value.elim <- as.numeric(ggdata_up$p.value.elim)


g1 <- ggplot(ggdata_up,
       aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(p.value.elim))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkred') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Heat stress responsive genes',
    subtitle = 'Upregulated genes - top 10 BP terms')+
  
  theme_bw(base_size = 22) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 10, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 10, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 10, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 15, face = "bold"), # Text size
    title = element_text(size = 15, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()




###downreuglated genes
ggdata_down <- ggdata[ggdata$regulation == 'down',]

ggdata_down <- ggdata_down %>%
  arrange(p.value.elim) %>%  # arrange in descending order
  slice(1:10) ##top10    

ggdata_down$Term <- factor(ggdata_down$Term, levels = rev(ggdata_down$Term)) # fixes order
ggdata_down$p.value.elim <- as.numeric(ggdata_down$p.value.elim)


g2 <- ggplot(ggdata_down, aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(p.value.elim))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkblue') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Heat stress responsive genes',
    subtitle = 'Downregulated genes - top 10 BP terms')+
  
  theme_bw(base_size = 22) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 10, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 10, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 10, face = 'bold'),
    axis.title.x = element_text(size = 10, face = 'bold'),
    axis.title.y = element_text(size = 10, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 15, face = "bold"), # Text size
    title = element_text(size = 15, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()

pdf(file = "/Users/mariebrasseur/Desktop/heatstress_upregulated_top10.pdf", width = 9, height = 6.5)
g1
dev.off()


pdf(file = "/Users/mariebrasseur/Desktop/heatstress_downregulated_top10.pdf", width = 10.5, height = 7)
g2
dev.off()


                     