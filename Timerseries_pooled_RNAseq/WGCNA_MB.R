###WGCNA
library("DESeq2")
library('sva')
#BioNERO package?
library('WGCNA')  
library('tximport')
library('readr')
library('PCAtools')
library('tidyr')
library('stringr')

##perform with data_adj (which is SVA corrected and variance stabilized)?? -> WGCNA at least suggets to remove batch effects
##accordign to https://support.bioconductor.org/p/42615/ a frozen sva wiould to the job

##setup
allowWGCNAThreads(30)
temp_cor <- cor       
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)
options(stringsAsFactors = FALSE)

###load data with tximport
#set to ballgown coverage data dir
setwd('../genome_guided_expression_gff_as_ref/')
files = c(read.table('samples_lst.txt')) ##file with .ctab file names
files_names <- c(files$V1)
files <- c(files$V2)
names(files) <- files_names

#create tx2gene file
a <- read_tsv(files[1]) ##somehow does not work when I want to store it in temp :/
tx2gene <- a[, c("t_name", "gene_id")]

#thisdoes not contain mitogenes
#tx2gene <- read.table('/share/pool/mbrasseur/eelpout_RNA_timeseries/reference/zoarces_tx2gene')

#read in tx abundances; for timeseries -> 150 bp
txi.stringtie <- tximport(files, type = "stringtie", tx2gene = tx2gene, readLength = 150)

#change to analysis dir
setwd('/share/pool/mbrasseur/eelpout_RNA_timeseries/WGCNA')

##load metadata
sampleData <- read.csv('../DESeq/timeseries_metadata.csv', sep = ';', header = 1)
row.names(sampleData) <- sampleData$Sample
sampleData <- sampleData[colnames(txi.stringtie$counts),]

sampleData$Year <- factor(sampleData$Year, levels = c(1994, 1999, 2005,2010,2012,2017,2019,2021)) ###treatment should be factor
levels(sampleData$Year) ##check the levels -> control must come first

sampleData$Location <- factor(sampleData$Location, levels = c("Dar", "Mel", "VaM"))

##check; if TRUE, proceed
all(colnames(txi.stringtie$counts) == rownames(sampleData)) 

#create DESeq object
dds1 <- DESeqDataSetFromTximport(txi.stringtie, colData=sampleData, design= ~ Location + Year) ##could make sense to fit time as numerical predictor an the add interaction term 

dds1 <- dds1[setdiff(rownames(dds1), '.'),]


##some filtering steps; omit this for now
dds1 <- estimateSizeFactors(dds1)
nc <- counts(dds1, normalized=TRUE)

filter2 <- rowSums(nc >= 3) >= 8 
dds1 <- dds1[filter2,]

##how many genes
length(row.names(counts(dds1))) ## 20041 genes


##run DESeq() on the complete data set and use the LRT
dds1 <- DESeq(dds1, test="LRT", reduced = ~ Location) ##this will test if all time factor levels significantly contribute to the likelihood of the mdodel

##correct for batch effect
dat  <- counts(dds1, normalized = TRUE)

##estimate latent factors i.e., surrogates
mod  <- model.matrix(~ Location + Year, colData(dds1)) ##specify full model
mod0 <- model.matrix(~   1, colData(dds1)) ##specify reduced model

##run sva
svseq <- svaseq(dat, mod, mod0) ##note down here the number of surrogates; 4


##frozen SVA to regress out latent factors in the PCA
newV = NULL ## neccessary due to bug in the sva pacakge
dat_vst  <- vst(dds1, blind=FALSE)

#run frozen sva
fsvaobj <- fsva(assay(dat_vst), mod, svseq, method = 'exact') 

##get the adjusted data
mat <- fsvaobj$db
head(mat)
##outlier genes and samples need to be removed -> different ways to identify these
##identify outlier genes with WGCAN function


##rows must be samples and cols must be genes
mat <- t(mat)

##check if there are too many missing values
gsg <- goodSamplesGenes(mat)
summary(gsg)
gsg$allOK ##if TRUE, then proceed

##if FALSE, then you need to filter out outliers
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(mat)[!gsg$goodGenes], collapse = ", "))); #Identifies and prints outlier genes
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(mat)[!gsg$goodSamples], collapse = ", "))); #Identifies and prints oulier samples
  mat <- mat[gsg$goodSamples == TRUE, gsg$goodGenes == TRUE] # Removes the offending genes and samples from the data
}

##evaluate again
gsg <- goodSamplesGenes(mat)
summary(gsg)
gsg$allOK 


##remove outlier samples based on hierarchical clustering
sampleTree <- hclust(dist(mat, method = 'can'), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
png('sample_clustering_outlier_detection.png')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()

##if outliers are identified, they might be removed with a cutree function 
#cutoff_value <- 72         ##this must be adjusted depending on the tree 

#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#     cex.axis = 1.5, cex.main = 2)
#draw on line to show cutoff height
#abline(h = cutoff_value, col = "red")


##cutoff
#cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = cutoff_value, minSize = 10) #returns numeric vector; Only branches whose size is at least minSize are retained
#Remove outlier
#mat <- mat[cut.sampleTree==1, ] we do not remove any outlier here because there are none


####network construction, 
###the first step is always to determine the power (?)

?pickSoftThreshold

spt <- pickSoftThreshold(mat, networkType = "signed", powerVector = seq(1, 40, by = 2), RsquaredCut = 0.85)

##plot R2 values as a function of the soft thresholds
par(mar=c(1,1,1,1))
png('r2_vs_power.png')
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.85,col="red")
dev.off()

##the best soft power threshold is the one which maximizes R2 (R^2 >=85 (or 80)??)
#and mean connectivity/minimizing the number of connections lost -> here 10

#Plot mean connectivity as a function of soft thresholds
par(mar=c(1,1,1,1))
png('mean_connectivity_vs_power.png')
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")
dev.off()


##NOTE: the higher the value, the stronger the connection strength 
#of highly correlated gene expression profiles -> the more devalued are low correlations

softPower <- 9



###the following code describes the MANUAL step-by-step network construction; adapted from: https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html

#1st step: identify some sort of similarity measurement between pairs of genes. -> similarity measurement represents the concordance of gene expression profiles across samples. 
#In WGCNA the similarity measurement for each pair of genes (gene i and gene j) is denoted by their Pearson correlation coefficient.
#for unsigned networks, you take the absolute value of the correlation, but signed networks are more common & biological interpretation is more straightforward


#2nd step: translate similarity measurement into the gene pairs adjacency to generate an adjacency matrix.
#Adjacency = assignment of a connection strength based on the co-expression similarity measurement (i.e. Pearson correlation coefficient). 
#Nodes are considered connected if they have a significant pairwise correlation association.

##adjacency calculation depends on the choice for unweighted vs weighted network 
#(unweighted networks are binary, whereas weighted networks can depict gradients in association
#weighted networks utilize a power function based on a soft threshold parameter ? which must be determined (determines the sensitivity and specificity of the analysis)


##create adjacency matrix for module construction 
?adjacency
adjacency <- adjacency(mat, power = softPower, type = "signed", distOptions = "method = 'canberra'")

##module construction;to assess topological overlap, adjacency must be transformed to dissimilarity;
#use TOM-based dissimilarity (topological overlap matrix)
?TOMsimilarity
TOM <- TOMsimilarity(adjacency, TOMType = "signed", verbose = 2) #this gives similarity
TOM.dissimilarity <- 1-TOM #this gives dissimilarity

###this dissimilarity can be used for gene clustering/finding modules = group of gene profiles that are highly correlated or have a high topological overlap.
##hierarchical clustering analysis

#The dissimilarity/distance measures are clustered and 
#a dendrogram (cluster tree) of genes is constructed.

geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
#plot dendrogram
sizeGrWindow(12,9)
png('genes_clustered_based_on_TOM.png')
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)
dev.off()


###identify modules
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, minClusterSize = 30, deepSplit=2)

#?cutreeDynamic ###params??

table(Modules) #returns a table of the counts of factor levels in an object.here: no. of genes assigned to each created module, whereby 0 is reserved for genes which do not fit into any module 

#plot module assignment under gene dendro
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors)
png('gene_dendro_with_assigned_modules.png')
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



###visualize the topological overlap of the network 
#apply a power for visualization 
plotTOM = TOM.dissimilarity^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function but reverse the colormapping
library(gplots)
#library(devEMF)

myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')

#options(expressions = 50000)

png('TOMplot.png', width = 7, height = 7, units = 'cm', res = 300)
TOMplot(plotTOM, geneTree, ModuleColors, main = "", col=myheatcol)
dev.off()


##identify eigengenes; 
MElist <- moduleEigengenes(mat, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)

#####module clustering based on eigengene similarity
#first: define distance between eigengenes
ME.dissimilarity = 1-cor(MEs)

METree = hclust(as.dist(ME.dissimilarity), method = "average") #Clustering eigengenes 
par(mar = c(0,4,2,0)) #seting margin sizes
par(cex = 0.6);#scaling the graphic
png('eigengen_clustering.png')
plot(METree)
abline(h=.2, col = "red") #a height of .25 corresponds to correlation of .75
dev.off()


##merge everything with < 0.2 dissimilarity
merge <- mergeCloseModules(mat, ModuleColors, cutHeight = .2)
# The merged module colors, assigning one color to each module
mergedColors = merge$colors
# Eigengenes of the new merged modules
mergedMEs = merge$newMEs

###compare merged vs unmerged modules
png("Dendrogram.png")
plotDendroAndColors(geneTree, cbind(ModuleColors, mergedColors), 
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules")
dev.off()


###plot how many genes per merged module
png("NumberOfGenesPerModule_merged.png")
barplot(table(mergedColors), las=2,cex.names = 0.8, main = "Number of genes per module after merging")
dev.off()

###plot how many genes per unmerged module
png("NumberOfGenesPerModule_unmerged.png")
barplot(table(ModuleColors), las=2,cex.names = 0.8, main = "Number of genes per module without merging")
dev.off()

###leave it unmerged since module merging does not really change anything


###perform correlation analyses
##first: biometric parameters i.e., average weight & body length

##check if the biometric data are normally distributed
setwd('C:/Users/mbras/Desktop/UPB - poolseq time line/')
biometrics <- read.csv('Eelpout_metadata_biometrics_and_time.csv', sep = ';')

colnames(biometrics)
##make Q-Q-plots
qqnorm(biometrics$Aalmutter.Gewicht.DAR)
qqline(biometrics$Aalmutter.Gewicht.DAR)
##test formally
shapiro.test(biometrics$Aalmutter.Gesamtl?nge.DAR)

##acc. to shapiro test, these two fail to be normally distributed
#Aalmutter.Gesamtl?nge.DAR
#Aalmutter.Gewicht.DAR



##plot with ggplot
library('ggplot2')
library('plyr')
library('gridExtra')

#transform to longformat

biometrics.long <- biometrics[,c(1:4, 8:14)]

biometrics.long <- gather(biometrics.long, biometric, value, EXACT_predicted_water_temp_at_sea_floor:Average_body_length)

class(biometrics.long$value)

unique(biometrics.long$biometric)

###biometrics

p0 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_body_length', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  #labs(title = 'Average body length')+
  ylab('Average body length [cm]')
ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_body_length', ]$value))

p0

p1 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_body_length', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = F) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Average body length [cm]')
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_body_length', ]$value))

p1 

p2 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_body_weight', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = F) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Average body weight [g]')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_body_weight', ]$value))

p2

p3 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_liver_weight', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = F) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Average liver weight [g]')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_liver_weight', ]$value))

p3


pdf("biometrics.pdf", width = 9, height = 7)
grid.arrange(p1, p2, p3, p0, ncol=2)
dev.off()

####env data

p4 <- ggplot(biometrics.long[biometrics.long$biometric == 'EXACT_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  labs(title = 'Predicted water temp. at sea floor')+
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Exact temperature [°C]')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'EXACT_predicted_water_temp_at_sea_floor', ]$value))

p4



p5 <- ggplot(biometrics.long[biometrics.long$biometric == 'MEDIAN_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value)) +
  geom_line(aes(color = Marine_region), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Median temperature [°C]')+
  labs(title = 'Predicted water temp. at sea floor')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'MEDIAN_predicted_water_temp_at_sea_floor', ]$value))

p5



p6 <- ggplot(biometrics.long[biometrics.long$biometric == 'MAX_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value)) +
  geom_line(aes(color = Marine_region), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Max. temperature [°C]')+
  labs(title = 'Predicted water temp. at sea floor')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'MAX_predicted_water_temp_at_sea_floor', ]$value))

p6


p7 <- ggplot(biometrics.long[biometrics.long$biometric == 'MIN_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value)) +
  geom_line(aes(color = Marine_region), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Min. temperature [°C]')+
  labs(title = 'Predicted water temp. at sea floor')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'MIN_predicted_water_temp_at_sea_floor', ]$value))

p7


pdf("env_data.pdf")
grid.arrange(p4, p5, p6, p7, ncol=2)
dev.off()




##load on cluster
allTraits <- read.csv("Eelpout_metadata_biometrics_and_time.csv", sep = ';')

##match trait data the expression data by sample ID, which must be also rownames
row.names(allTraits) <- allTraits$ID
allTraits <- allTraits[row.names(mat),]

colnames(allTraits)

##subset to env data
envTraits <- allTraits[, c(1, 8:14)]
envTraits$Year <- envTraits$Year - 1994

# Define numbers of genes and samples
nGenes = ncol(mat)
nSamples = nrow(mat) ##this must be adjusted for biometric data


##calculate eigengene gene significance
##if we want to use unmerged MEs, then use the object MEs

#make sure that the right cor function is used here! 
module.trait.correlation.env = cor(MEs, envTraits, method = "pearson") 
module.trait.Pvalue.env = corPvalueStudent(module.trait.correlation.env, nSamples) #calculate the p-value associated with the correlation


module.trait.Pvalue.env.adj <- module.trait.Pvalue.env

'''diese reihenfolge:
[1] "Year"                                    
[2] "EXACT_predicted_water_temp_at_sea_floor" 
[3] "MEDIAN_predicted_water_temp_at_sea_floor"
[4] "MAX_predicted_water_temp_at_sea_floor"   
[5] "MIN_predicted_water_temp_at_sea_floor"   
[6] "Average_body_weight"                     
[7] "Average_liver_weight"                    
[8] "Average_body_length" '''



module.trait.Pvalue.env.adj[,1] <- p.adjust(module.trait.Pvalue.env, method = "BH")[1:23] #23 MEs
module.trait.Pvalue.env.adj[,2] <- p.adjust(module.trait.Pvalue.env, method = "BH")[24:46]
module.trait.Pvalue.env.adj[,3] <- p.adjust(module.trait.Pvalue.env, method = "BH")[47:69]
module.trait.Pvalue.env.adj[,4] <- p.adjust(module.trait.Pvalue.env, method = "BH")[70:92]
module.trait.Pvalue.env.adj[,5] <- p.adjust(module.trait.Pvalue.env, method = "BH")[93:115]
module.trait.Pvalue.env.adj[,6] <- p.adjust(module.trait.Pvalue.env, method = "BH")[116:138]
module.trait.Pvalue.env.adj[,7] <- p.adjust(module.trait.Pvalue.env, method = "BH")[139:161]
module.trait.Pvalue.env.adj[,8] <- p.adjust(module.trait.Pvalue.env, method = "BH")[162:184]



module.trait.Pvalue.env.adj <- as.data.frame(module.trait.Pvalue.env.adj)
module.trait.correlation.env <- as.data.frame(module.trait.correlation.env)

module.trait.correlation.env$Year[module.trait.Pvalue.env.adj$Year > 0.05] <- NA
module.trait.correlation.env$EXACT_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$EXACT_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.correlation.env$MEDIAN_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$MEDIAN_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.correlation.env$MAX_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$MAX_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.correlation.env$MIN_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$MIN_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.correlation.env$Average_body_weight[module.trait.Pvalue.env.adj$Average_body_weight > 0.05] <- NA
module.trait.correlation.env$Average_liver_weight[module.trait.Pvalue.env.adj$Average_liver_weight > 0.05] <- NA
module.trait.correlation.env$Average_body_length[module.trait.Pvalue.env.adj$Average_body_length > 0.05] <- NA

#module.trait.correlation.env <- na.omit(module.trait.correlation.env)


module.trait.Pvalue.env.adj$Year[module.trait.Pvalue.env.adj$Year > 0.05] <- NA
module.trait.Pvalue.env.adj$EXACT_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$EXACT_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.Pvalue.env.adj$MEDIAN_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$MEDIAN_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.Pvalue.env.adj$MAX_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$MAX_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.Pvalue.env.adj$MIN_predicted_water_temp_at_sea_floor[module.trait.Pvalue.env.adj$MIN_predicted_water_temp_at_sea_floor > 0.05] <- NA
module.trait.Pvalue.env.adj$Average_body_weight[module.trait.Pvalue.env.adj$Average_body_weight > 0.05] <- NA
module.trait.Pvalue.env.adj$Average_liver_weight[module.trait.Pvalue.env.adj$Average_liver_weight > 0.05] <- NA
module.trait.Pvalue.env.adj$Average_body_length[module.trait.Pvalue.env.adj$Average_body_length > 0.05] <- NA

#module.trait.Pvalue.env.adj <- na.omit(module.trait.Pvalue.env.adj)



# Will display correlations and their p-values; input must be matrix not dataframe
textMatrix = paste(signif(as.matrix(module.trait.correlation.env), 2), "\n(",
                   signif(as.matrix(module.trait.Pvalue.env.adj), 1), ")", sep = "")

dim(textMatrix) = dim(module.trait.correlation.env)


# Display the correlation values within a heatmap plot
pdf('module_trait_pearson_correlation_adjusted_p.pdf', height = 7, width = 7)
labeledHeatmap(Matrix = module.trait.correlation.env,
               xLabels = c("Year", "Exact temp.","Median seasonal temp.", "Max. seasonal temp.", "Min. seasonal temp.", "Average body weight", "Average liverweight", "Average body length"),
               yLabels = row.names(module.trait.correlation.env),
               ySymbols = row.names(module.trait.correlation.env),
               #yLabels = names(MEs),
               #ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 0.4,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


table(ModuleColors)

cor <- temp_cor     # Return cor function to original namespace



#####plot expression profiles & eigengenes

lvls =  c( "Dar_94", "Dar_99", "Dar_05", "Dar_10", "Dar_12", "Dar_17", "Dar_19", "Dar_21",
           "Mel_94", "Mel_99", "Mel_05", "Mel_10", "Mel_12", "Mel_17", "Mel_19", "Mel_21", 
           "VaM_94", "VaM_99", "VaM_05", "VaM_10", "VaM_12", "VaM_17", "VaM_19", "VaM_21")


###all MEs 
ME.plot <- MEs
ME.plot$name <- factor(row.names(ME.plot), levels = lvls)


##transform to long format with tidyr
ME.plot.long <- gather(ME.plot, module, expression, MEblack:MEyellow)

pdf('MEs_expression_profiles_all.pdf')
ggplot(ME.plot.long, aes(x=name, y=expression, group=module)) +
  geom_line(aes(color = module), show.legend = FALSE) +
  scale_color_manual(values = str_remove(colnames(ME.plot), 'ME')) +
  theme_light() +
  theme(text = element_text(size = 7.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  facet_wrap(~ module, ncol= 3) +
  labs(x = NULL,
       y = "Eigengene expression")
dev.off()



###all MEs that had a sig. correlation to any of the tested traits
ME.sub <- MEs[, -c(9,10,18,20)] ##removal of grey, grey60, red, salmon
ME.sub$name <- factor(row.names(ME.sub), levels = lvls)


##transform to long format with tidyr
ME.sub.long <- gather(ME.sub, module, expression, MEblack:MEyellow)

pdf('MEs_expression_profiles_all_with_corr.pdf')
ggplot(ME.sub.long, aes(x=name, y=expression, group=module)) +
  geom_line(aes(color = module), show.legend = FALSE) +
  scale_color_manual(values = str_remove(colnames(ME.sub), 'ME')) +
  theme_light() +
  theme(text = element_text(size = 7.5), axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)) +
  facet_wrap(~ module, ncol= 3) +
  labs(x = NULL,
       y = "Eigengene expression")
dev.off()



###pick out all MEs with significant env correlation 

modules.env <- row.names(module.trait.correlation.env)

ME.env <- as.data.frame(MEs[, modules.env])
ME.env$name <- row.names(ME.env)
ME.env$name <- factor(ME.env$name, levels = lvls)


##transform to long format with tidyr
ME.env.long <- gather(ME.env, module, expression, MEblue:MEyellow)

pdf('MEs_env_expression_profiles.pdf')
ggplot(ME.env.long, aes(x=name, y=expression, group=module)) +
  geom_line(aes(color = module), show.legend = FALSE) +
  scale_color_manual(values = str_remove(modules.env, 'ME')) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Sample",
       y = "Eigengene expression")
dev.off()


###biometric modules
modules.biom <- row.names(module.trait.correlation.biom)

ME.biom <- as.data.frame(MEs[, modules.biom])

ME.biom$name <- row.names(ME.biom)
ME.biom$name <- factor(ME.biom$name, levels = lvls)


##transform to long format
ME.biom.long <- gather(ME.biom, module, expression, MEblack:MEturquoise)

pdf('MEs_biom_expression_profiles.pdf')
ggplot(ME.biom.long, aes(x=name, y=expression, group=module)) +
  geom_line(aes(color = module), show.legend = FALSE) +
  scale_color_manual(values = str_remove(modules.biom, 'ME')) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Sample",
       y = "Eigengene expression")
dev.off()




#####################################
#individual genes
#modGenes = (ModuleColors==module)
modGenes=mat[,ModuleColors==module] ##this is the normalized expression of the genes in the blue module
#dimensions are samples x genes


module_df <- data.frame(
  gene_id = colnames(modGenes),
  colors = module
)


# Pull out list of genes in that module
library(dplyr)

submod_df = data.frame(t(modGenes)) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    colors = module
  )

##check if this is right
png('test.png')
ggplot(submod_df, aes(x=factor(name, levels = lvls), y=scale(value, scale = TRUE), group=gene_id)) +
  geom_line(aes(color = colors),
            alpha = 0.05) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
#  facet_grid(rows = vars(module)) +
  labs(x = "sample",
       y = "normalized expression")
dev.off()



################################


##plot the module MEs which are significantly correlated with both, biometric data and temperature vars

##all in env except yellow
modules.both <- row.names(module.trait.correlation.env)[1:5]


ME.both <- as.data.frame(MEs[, modules.both])
ME.both$name <- row.names(ME.both)
ME.both$name <- factor(ME.both$name, levels = lvls)


##transform to long format with tidyr
ME.both.long <- gather(ME.both, module, expression, MEblue:MEturquoise)

colors.both <- gsub('MEgreen', 'MEdarkgreen', modules.both)
colors.both <- str_remove(colors.both, 'ME')

module_names <- c(
  `MEblue` = "Blue module",
  `MEgreen` = "Green module",
  `MElightgreen` = "Lightgreen module",
  `MEmidnightblue` = "Midnightblue module",
  `MEturquoise` = "Turquoise module"
)


pdf('MEs_both_vars_expression_profiles.pdf')
ggplot(ME.both.long, aes(x=name, y=expression, group=module)) +
  geom_line(aes(color = module), show.legend = FALSE) +
  scale_color_manual(values = colors.both) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module),labeller = as_labeller(module_names)) +
  labs(x = "Sample",
       y = "Eigengene expression")
dev.off()


####plot correlation of eigengene expression and traits

##env data is in envTraits
##biometric data is in biometricTraits

#row.names(ME.both) == row.names(envTraits) TRUE
#row.names(ME.both) == row.names(biometricTraits) TRUE

##select only one  trait with the highest correlation on average;

#min. annual temp.
ME.both$min_temp <- envTraits$MIN_predicted_water_temp_at_sea_floor

#average body length
ME.both$length <- biometricTraits$Average_body_length

ME.both.long <- gather(ME.both, module, expression, MEblue:MEturquoise)

pdf('Corr_MEexpression_min_temp.pdf')
ggplot(ME.both.long, aes(x=min_temp, y=expression, label= name)) +
  geom_point(aes(color = module), show.legend = FALSE) +
  geom_text(hjust=0, vjust=0) +
  scale_color_manual(values =colors.both) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Min. temperature",
       y = "Eigengene expression")
dev.off()


pdf('Corr_MEexpression_body_length.pdf')
ggplot(ME.both.long, aes(x=length, y=expression, label= name)) +
  geom_point(aes(color = module), show.legend = FALSE) +
  geom_text(hjust=0, vjust=0) +
  scale_color_manual(values =colors.both) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "Average body length",
       y = "Eigengene expression")
dev.off()

###plot temperature and body length

png('body_length.png')
ggplot(ME.both.long, aes(x=name, y=length)) +
  geom_point()+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75))
dev.off()


png('min_temp.png')
ggplot(ME.both.long, aes(x=name, y=min_temp)) +
  geom_point()+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75))
dev.off()




##find overlaps between the different data sets

##create all the gene lists

heatstress <- read.csv('../../eelpout_RNA_heatstress_experiment/DESeq/DEG_zoarces_heatstress_exp.csv', row.names = 1)

heatstress$padj[is.na(heatstress$padj)] <- 1
heatstress <- heatstress[heatstress$padj < 0.05,]
heatstress <- row.names(heatstress)


time.series <- read.csv('../DESeq/DEG_zoarces_timeseries_LRT.csv', row.names = 1)

time.series$padj[is.na(time.series$padj)] <- 1
time.series <- time.series[time.series$padj < 0.05,]
time.series <- row.names(time.series)

##modules
allModules <- names(table(ModuleColors))


###first; check that no genes is (erroneously) in two modules
for (col in allModules){
  for (i in c(1:length(allModules))){
    geneset1 = colnames(mat[,ModuleColors==col])
    geneset2 = colnames(mat[,ModuleColors==allModules[i]])
    intersection = intersect(geneset1, geneset2)
    if (length(intersection > 0)){
      print(paste('shared genes between ', toString(col), ' and ', toString(allModules[i]), ' : ', toString(length(intersection))))
    }
  }
}

######
## initialize the dataframe; required format is:
## | FROM | TO | VALUE

df = data.frame(
  FROM = c('heatstress'),
  TO = c('timeseries'),
  VALUE = length(intersect(heatstress, time.series)))


for (col in allModules){
    #find intersection between heatstress and modules
    geneset_temp = colnames(mat[,ModuleColors==col])
    intersection = intersect(heatstress, geneset_temp)
    new_row = list(FROM = 'heatstress', TO = col, VALUE = length(intersection))
    df = rbind(df, new_row)
    
    #find intersection between timeseries and modules
    intersection = intersect(time.series, geneset_temp)
    new_row = list(FROM='timeseries', TO= col, VALUE = length(intersection))
    df = rbind(df, new_row)
}




#####now find the genes which are unique
all = c()

for (col in allModules){
  all = c(all, colnames(mat[,ModuleColors==col]))
}

###20041 genes in these modules

all.time = unique(c(all, time.series)) ##still 20041

setdiff(heatstress, all.time) #one gene is unique to heatstressm and not in all others

all.heatstress = unique(c(all, heatstress)) #20042

setdiff(time.series, all.heatstress) #time series completely included in all others

df = rbind(df, list(FROM = 'heatstress', TO = 'heatstress', VALUE = length(setdiff(heatstress, all.time))))

###unique ones from the modules!

all.deseq = unique(c(heatstress, time.series))

for (col in allModules){
  geneset_temp = colnames(mat[,ModuleColors==col])
  unique = setdiff(geneset_temp, all.deseq)
  new_row = list(FROM = col, TO = col, VALUE = length(unique))
  df = rbind(df, new_row)
}



###this df can now be used for circlize chord diagram

library(circlize)

length(allModules) #23

group <- structure(c(rep('DESeq2', each = 2),rep('WGCNA', each=23)), names =  c('heatstress', 'timeseries', allModules))

#grid.col = c('blue', 'red',rep('grey', each=23))
grid.col = c('#6F4D38', '#C2A77C', allModules)


##for link color mapping
border <- rep(0, each = length(row.names(df)))

idx.heatstress <- which(df$FROM == 'heatstress') 
idx.timeseries <- which(df$FROM == 'timeseries')

border[idx.heatstress] <- '#6F4D38'
border[idx.timeseries] <- '#C2A77C'

circos.par(start.degree = -77)
pdf('circos_plot_intersection_and_unique_genes_scaled.pdf')
chordDiagram(df,group = group, grid.col= grid.col, col='grey', transparency = 0.5, 
             link.border = border,
             link.zindex = rev(c(1:length(row.names(df)))),
             #self.link = 2, 
             scale = TRUE,
             annotationTrack = "grid", big.gap=15)
dev.off()
circos.clear()




############now plot only the intersections
df = data.frame(
  FROM = c('heatstress'),
  TO = c('timeseries'),
  VALUE = length(intersect(heatstress, time.series)))


for (col in allModules){
  #find intersection between heatstress and modules
  geneset_temp = colnames(mat[,ModuleColors==col])
  intersection = intersect(heatstress, geneset_temp)
  new_row = list(FROM = 'heatstress', TO = col, VALUE = length(intersection))
  df = rbind(df, new_row)
  
  #find intersection between timeseries and modules
  intersection = intersect(time.series, geneset_temp)
  new_row = list(FROM='timeseries', TO= col, VALUE = length(intersection))
  df = rbind(df, new_row)
}


grid.col = c('#6F4D38', '#C2A77C', allModules)
#grid.col = c('#C2A77C', '#C2A77C', allModules)

col_fun = colorRamp2(range(df$VALUE), c("#C0C0C0", "#708090"), transparency = 0.5)

circos.par(start.degree = -12)
pdf('circos_plot_only_intersections.pdf')
chordDiagram(df,group = group, grid.col= grid.col, col=col_fun, transparency = 0.5, 
             #link.border = c(rep('orange', each= 48), rep(0, each= 23)),
             #link.zindex = rank(df[[3]]),
             #self.link = 2, scale = TRUE,
             annotationTrack = "grid", 
             big.gap=15)
dev.off()
circos.clear()



###functional enrichment 
'''
test the genes in the following modules against the whole timeseries transcriptome:
blue
green
lightgreen
midnightblue
turquoise
yellow
'''


#################functional enrichment:
library(tidyr)
library(GO.db)
library(topGO)

setwd('../functional annotation')

#######if gene2go is already prepared
##check that no invisible line breaks are in the file ;open the file one time a just save it
#transform to longformat
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

background <- colnames(mat)

'''
blue
green
lightgreen
midnightblue
turquoise
yellow
'''

col='blue'


geneList <- colnames(mat[,ModuleColors==col]) ##this must get the ones
NOTgene_list <- setdiff(background, geneList) ##this must get the zeros

temp <- c(rep(1, each = length(geneList)), rep(0, each = length(NOTgene_list)))
names(temp) <- c(geneList, NOTgene_list)

##Create topGOdata object:

GOdata <- new("topGOdata",
                   ontology = "BP",  #ontology criteria = biological process
                   allGenes = temp, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                   geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

##run enrichment test: here Fishers Exact Test
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultFisher.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

#summarize in table
GOtable <- GenTable(GOdata, p.value.elim = resultFisher.elim, p.value.weight01 = resultFisher.weight01, topNodes = length(resultFisher.weight01@score), numChar = 10000000)

write.csv(GOtable, file = paste(col, '_WGCNA_module_BP_enrichment_fisher.csv', sep = ''))

#write.csv(GOtable, file = paste(col, '_WGCNA_module_MF_enrichment_fisher.csv', sep = ''))













############visualization

##module blue

ggdata <- read.csv("/Users/mariebrasseur/Desktop/trier/WGCNA/functional enrichment/blue_WGCNA_module_BP_enrichment_fisher.csv", sep = ';', row.names = 1)
ggdata <- down[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata$weight01 <- as.numeric(ggdata$weight01)

pdf(file= 'GOenrichment_pesticide_under_bi_downregulated_BP_top15.pdf', height = 9, width = 14)
gg1 <- ggplot(ggdata, aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(weight01))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkblue') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Insecticide effect L. basale -\ndownregulated genes under biotic interaction',
    subtitle = 'Top 15 Biological Process terms')+
  
  theme_bw(base_size = 28) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 19, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 15, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 15, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 15, face = 'bold'),
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), #removes the border
    legend.key.size = unit(1, "cm"), #Sets overall area/size of the legend
    legend.text = element_text(size = 19, face = "bold"), # Text size
    title = element_text(size = 19, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()
gg1
dev.off()


#upregulated genes
ggdata <- up[1:ntop,]
ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
ggdata$weight01 <- as.numeric(ggdata$weight01)

pdf(file= 'GOenrichment_pesticide_under_bi_upregulated_BP_top15.pdf',  height = 9, width = 14)
gg1 <- ggplot(ggdata,
              aes(x = Term, y = Significant, size = Significant/Expected, fill = -log10(weight01))) +
  expand_limits(y = 7) +
  geom_point(shape = 21) +  #shape = 21 allows to fill with other color
  scale_size(range = c(1,10)) +#, breaks = c(1, 4, 10)) +
  scale_fill_continuous(low = 'yellow', high = 'darkred') +
  guides(size=guide_legend("Enrichment ratio"))+ 
  
  xlab('') + ylab('Annotated genes') + 
  labs(
    title = 'Insecticide effect L. basale -\nupregulated genes under biotic interaction',
    subtitle = 'Top 15 Biological Process terms')+
  
  theme_bw(base_size = 28) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 19, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 15, face = 'bold', vjust = 1),
    
    axis.text.x = element_text(angle = 0, size = 15, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 15, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 15, face = 'bold'),
    axis.title.x = element_text(size = 15, face = 'bold'),
    axis.title.y = element_text(size = 15, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    
    #Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 19, face = "bold"), # Text size
    title = element_text(size = 19, face = "bold")) +
  
  #change legend title
  labs(fill = '-log10(p-value)') +
  
  coord_flip()
gg1
dev.off()
