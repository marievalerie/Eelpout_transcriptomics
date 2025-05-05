library("DESeq2")
library('sva')
library('WGCNA')  
library('tximport')
library('readr')
library('PCAtools')
library('tidyr')
library('stringr')



###########################
#WGCNA
###########################


#perform the WGCNA based on dcount data that is batch effect adjusted, if one exist i.e., SVA corrected and variance stabilized counts
#accordign to https://support.bioconductor.org/p/42615/ a frozen sva wiould to the job

#setup
allowWGCNAThreads(30)
temp_cor <- cor       
cor <- WGCNA::cor # Force it to use WGCNA cor function (fix a namespace conflict issue)
options(stringsAsFactors = FALSE)

#load data with tximport
#set to ballgown coverage data dir
setwd('../genome_guided_expression_gff_as_ref/')
files = c(read.table('samples_lst.txt')) ##file with .ctab file names
files_names <- c(files$V1)
files <- c(files$V2)
names(files) <- files_names

#create tx2gene file
a <- read_tsv(files[1]) 
tx2gene <- a[, c("t_name", "gene_id")]
tx2gene <- tx2gene[tx2gene$gene_id != ".",] 

#read in tx abundances; for timeseries -> 150 bp
txi.stringtie <- tximport(files, type = "stringtie", tx2gene = tx2gene, readLength = 150)

#change to analysis dir
setwd('/share/pool/mbrasseur/eelpout_RNA_timeseries/WGCNA')

#load metadata
sampleData <- read.csv('../DESeq/timeseries_metadata.csv', sep = ';', header = 1)
row.names(sampleData) <- sampleData$Sample
sampleData <- sampleData[colnames(txi.stringtie$counts),]

sampleData$Year <- factor(sampleData$Year, levels = c(1994, 1999, 2005,2010,2012,2017,2019,2021)) ###treatment should be factor
levels(sampleData$Year) ##check the levels -> control must come first

sampleData$Location <- factor(sampleData$Location, levels = c("Dar", "Mel", "VaM"))

#check; if TRUE, proceed
all(colnames(txi.stringtie$counts) == rownames(sampleData)) 

#create DESeq object
dds1 <- DESeqDataSetFromTximport(txi.stringtie, colData=sampleData, design= ~ Location + Year) 

#remove mitochondrial tRNAs and rRNAs
dds1 <- dds1[setdiff(rownames(dds1), '.'),]

#some filtering steps & DESeq normalization
dds1 <- estimateSizeFactors(dds1)
nc <- counts(dds1, normalized=TRUE)

filter2 <- rowSums(nc >= 3) >= 8 
dds1 <- dds1[filter2,]
dds1 <- DESeq(dds1, test="LRT", reduced = ~ Location) 

####correct for batch effect
dat  <- counts(dds1, normalized = TRUE)

#estimate latent factors i.e., surrogates
mod  <- model.matrix(~ Location + Year, colData(dds1)) ##specify full model
mod0 <- model.matrix(~   1, colData(dds1)) ##specify reduced model

#run sva
svseq <- svaseq(dat, mod, mod0) ##note down here the number of surrogates; 4

#frozen SVA to regress out latent factors
newV = NULL ## neccessary due to bug in the sva pacakge
dat_vst  <- vst(dds1, blind=FALSE)

#run frozen sva
fsvaobj <- fsva(assay(dat_vst), mod, svseq, method = 'exact') 

#get the adjusted data as matrix and use this in WGCNA
mat <- fsvaobj$db
head(mat)


####WGCNA
#outlier genes and samples need to be removed -> different ways to identify these
#identify outlier genes with WGCNA function

#rows must be samples and cols must be genes
mat <- t(mat)

#check if there are too many missing values
gsg <- goodSamplesGenes(mat)
summary(gsg)
gsg$allOK #if TRUE, proceed


#remove outlier samples based on hierarchical clustering
sampleTree <- hclust(dist(mat, method = 'can'), method = "average") #Clustering samples based on distance 

#Setting the graphical parameters
par(cex = 0.6);
par(mar = c(0,4,2,0))

#Plotting the cluster dendrogram
png('sample_clustering_outlier_detection.png')
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


####Network construction, 
#1. determine the power

spt <- pickSoftThreshold(mat, networkType = "signed", powerVector = seq(1, 40, by = 2), RsquaredCut = 0.85)

#plot R2 values as a function of the soft thresholds
par(mar=c(1,1,1,1))
png('r2_vs_power.png')
plot(spt$fitIndices[,1],spt$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(spt$fitIndices[,1],spt$fitIndices[,2],col="red")
abline(h=0.85,col="red")
dev.off()

#the best soft power threshold is the one which maximizes R2 
#and minimizes the number of connections lost 

#Plot mean connectivity as a function of soft thresholds
par(mar=c(1,1,1,1))
png('mean_connectivity_vs_power.png')
plot(spt$fitIndices[,1], spt$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(spt$fitIndices[,1], spt$fitIndices[,5], labels= spt$fitIndices[,1],col="red")
dev.off()


#NOTE: the higher the value, the stronger the connection strength 
#of highly correlated gene expression profiles -> the more devalued are low correlations

softPower <- 9

####create adjacency matrix for module construction 
adjacency <- adjacency(mat, power = softPower, type = "signed", distOptions = "method = 'canberra'")

###module construction;to assess topological overlap, adjacency must be transformed to dissimilarity;
#use TOM-based dissimilarity (topological overlap matrix)

TOM <- TOMsimilarity(adjacency, TOMType = "signed", verbose = 2) #this gives similarity
TOM.dissimilarity <- 1-TOM #this gives dissimilarity

#clustering
geneTree <- hclust(as.dist(TOM.dissimilarity), method = "average") 
plot dendrogram
sizeGrWindow(12,9)
png('genes_clustered_based_on_TOM.png')
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04)
dev.off()


#identify modules
Modules <- cutreeDynamic(dendro = geneTree, distM = TOM.dissimilarity, minClusterSize = 30, deepSplit=2)
table(Modules) 

#plot module assignment under gene dendro
ModuleColors <- labels2colors(Modules) #assigns each module number a color
table(ModuleColors)
png('gene_dendro_with_assigned_modules.png')
plotDendroAndColors(geneTree, ModuleColors,"Module",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



####visualize the topological overlap of the network 
#apply a power for visualization 
plotTOM = TOM.dissimilarity^7

#set diagonal to NA for a nicer plot
diag(plotTOM) = NA

#call the plot function but reverse the colormapping
library(gplots)
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')

png('TOMplot.png', width = 7, height = 7, units = 'cm', res = 300)
TOMplot(plotTOM, geneTree, ModuleColors, main = "", col=myheatcol)
dev.off()


####Identify eigengenes; 
MElist <- moduleEigengenes(mat, colors = ModuleColors) 
MEs <- MElist$eigengenes 
head(MEs)


####correlate Eigengenes with external traits:
#biometric parameters (average weight & body length) & environmental data

#check if the biometric data are normally distributed
setwd('C:/Users/mbras/Desktop/UPB - poolseq time line/')
biometrics <- read.csv('Eelpout_metadata_biometrics_and_time.csv', sep = ';')

#make Q-Q-plots
qqnorm(biometrics$Aalmutter.Gewicht.DAR)
qqline(biometrics$Aalmutter.Gewicht.DAR)
#test formally
shapiro.test(biometrics$Aalmutter.Gesamtlänge.DAR)

#plot with ggplot
library('ggplot2')
library('plyr')
library('gridExtra')

#transform to longformat
biometrics.long <- biometrics[,c(1:4, 8:14)]

biometrics.long <- gather(biometrics.long, biometric, value, EXACT_predicted_water_temp_at_sea_floor:Average_body_length)

class(biometrics.long$value)
unique(biometrics.long$biometric)


p0 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_body_length', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  #labs(title = 'Average body length')+
  ylab('Average body length [cm]')
ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_body_length', ]$value))


p1 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_body_length', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = F) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Average body length [cm]')
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_body_length', ]$value))


p2 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_body_weight', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = F) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Average body weight [g]')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_body_weight', ]$value))

p3 <- ggplot(biometrics.long[biometrics.long$biometric == 'Average_liver_weight', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = F) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Average liver weight [g]')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'Average_liver_weight', ]$value))

pdf("biometrics.pdf", width = 9, height = 7)
grid.arrange(p1, p2, p3, p0, ncol=2)
dev.off()


#environmental data

p4 <- ggplot(biometrics.long[biometrics.long$biometric == 'EXACT_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value )) +
  geom_line(aes(color = Sampling_site), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkgreen', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  labs(title = 'Predicted water temp. at sea floor')+
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Exact temperature [°C]')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'EXACT_predicted_water_temp_at_sea_floor', ]$value))


p5 <- ggplot(biometrics.long[biometrics.long$biometric == 'MEDIAN_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value)) +
  geom_line(aes(color = Marine_region), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Median temperature [°C]')+
  labs(title = 'Predicted water temp. at sea floor')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'MEDIAN_predicted_water_temp_at_sea_floor', ]$value))


p6 <- ggplot(biometrics.long[biometrics.long$biometric == 'MAX_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value)) +
  geom_line(aes(color = Marine_region), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Max. temperature [°C]')+
  labs(title = 'Predicted water temp. at sea floor')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'MAX_predicted_water_temp_at_sea_floor', ]$value))


p7 <- ggplot(biometrics.long[biometrics.long$biometric == 'MIN_predicted_water_temp_at_sea_floor', ], aes(x=Year, y=value)) +
  geom_line(aes(color = Marine_region), show.legend = T) +
  scale_color_manual(values = c('purple', 'darkorange')) +
  scale_x_continuous('Year', labels = unique(biometrics$Year), breaks = unique(biometrics$Year))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1)) +
  ylab('Min. temperature [°C]')+
  labs(title = 'Predicted water temp. at sea floor')+
  ylim(0,max(biometrics.long[biometrics.long$biometric == 'MIN_predicted_water_temp_at_sea_floor', ]$value))


pdf("env_data.pdf")
grid.arrange(p4, p5, p6, p7, ncol=2)
dev.off()



#####correlate the traits
#load on cluster
allTraits <- read.csv("Eelpout_metadata_biometrics_and_time.csv", sep = ';')

#match trait data the expression data by sample ID, which must be also rownames
row.names(allTraits) <- allTraits$ID
allTraits <- allTraits[row.names(mat),]

#subset to env data
envTraits <- allTraits[, c(1, 8:14)]
envTraits$Year <- envTraits$Year - 1994

#define numbers of genes and samples
nGenes = ncol(mat)
nSamples = nrow(mat) ##this must be adjusted for biometric data

#make sure that the right cor function is used here! 
module.trait.correlation.env = cor(MEs, envTraits, method = "pearson") 
module.trait.Pvalue.env = corPvalueStudent(module.trait.correlation.env, nSamples) #calculate the p-value associated with the correlation

module.trait.Pvalue.env.adj <- module.trait.Pvalue.env

'''
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

#display correlations and their p-values; input must be matrix (no df)
textMatrix = paste(signif(as.matrix(module.trait.correlation.env), 2), "\n(",
                   signif(as.matrix(module.trait.Pvalue.env.adj), 1), ")", sep = "")

dim(textMatrix) = dim(module.trait.correlation.env)


#display the correlation values within a heatmap plot
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

cor <- temp_cor     # Return cor function to original namespace


####plot eigengenes (i.e., average expression profiles of modules) 

lvls =  c( "Dar_94", "Dar_99", "Dar_05", "Dar_10", "Dar_12", "Dar_17", "Dar_19", "Dar_21",
           "Mel_94", "Mel_99", "Mel_05", "Mel_10", "Mel_12", "Mel_17", "Mel_19", "Mel_21", 
           "VaM_94", "VaM_99", "VaM_05", "VaM_10", "VaM_12", "VaM_17", "VaM_19", "VaM_21")

ME.plot <- MEs
ME.plot$name <- factor(row.names(ME.plot), levels = lvls)

#transform to long format with tidyr
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


####find overlaps between the different data sets with circlize chord diagram

#create the gene lists

#1. heat stress responsive genes
heatstress <- read.csv('../../eelpout_RNA_heatstress_experiment/DESeq/DEG_zoarces_heatstress_exp.csv', row.names = 1)
heatstress$padj[is.na(heatstress$padj)] <- 1
heatstress <- heatstress[heatstress$padj < 0.05,]
heatstress <- row.names(heatstress)

#2. time responsive genes
time.series <- read.csv('../DESeq/DEG_zoarces_timeseries_LRT.csv', row.names = 1)
time.series$padj[is.na(time.series$padj)] <- 1
time.series <- time.series[time.series$padj < 0.05,]
time.series <- row.names(time.series)

#3. genes in WGCNA modules
allModules <- names(table(ModuleColors))

#just to be sure: check that no genes is (erroneously) in two modules
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

####merge the gene lists
#initialize the dataframe; required format is:
#| FROM | TO | VALUE

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


####nfind the genes which are unique
all = c()
for (col in allModules){
  all = c(all, colnames(mat[,ModuleColors==col]))
}


####20041 genes in the modules
all.time = unique(c(all, time.series)) #still 20041
setdiff(heatstress, all.time) #one gene is unique to heatstressm and not in all others
all.heatstress = unique(c(all, heatstress)) #20042
setdiff(time.series, all.heatstress) #time series completely included in all others

df = rbind(df, list(FROM = 'heatstress', TO = 'heatstress', VALUE = length(setdiff(heatstress, all.time))))

####unique genes from the WGCNA modules
all.deseq = unique(c(heatstress, time.series))

for (col in allModules){
  geneset_temp = colnames(mat[,ModuleColors==col])
  unique = setdiff(geneset_temp, all.deseq)
  new_row = list(FROM = col, TO = col, VALUE = length(unique))
  df = rbind(df, new_row)
}

#this df can now be used for circlize chord diagram

library(circlize)
length(allModules) #23

group <- structure(c(rep('DESeq2', each = 2),rep('WGCNA', each=23)), names =  c('heatstress', 'timeseries', allModules))

grid.col = c('#6F4D38', '#C2A77C', allModules)

#for link color mapping
border <- rep(0, each = length(row.names(df)))

idx.heatstress <- which(df$FROM == 'heatstress') 
idx.timeseries <- which(df$FROM == 'timeseries')

border[idx.heatstress] <- '#6F4D38'
border[idx.timeseries] <- '#C2A77C'

circos.par(start.degree = -77)
pdf('circos_plot_intersection_and_unique_genes.pdf')
chordDiagram(df,group = group, grid.col= grid.col, col='grey', transparency = 0.5, 
             link.border = border,
             link.zindex = rev(c(1:length(row.names(df)))),
             #self.link = 2, 
             #scale = TRUE,
             annotationTrack = "grid", big.gap=15)
dev.off()
circos.clear()



###########################
#Functional enrichment
###########################

library(tidyr)
library(GO.db)
library(topGO)


####test the genes in the two main modules blue and turquoise against the whole timeseries transcriptome:

setwd('../functional annotation')

#if gene2go is already prepared
#check that no invisible line breaks are in the file

GO_ids = read.csv('StringtieEggNog_Gene2GO.csv', sep=';', header = F)

#transform to longformat
long_GO <- gather(GO_ids, Gen, IDs, V2:V1160)

#take out genes without GO terms
long_GO <- long_GO[which(long_GO$`IDs` != ""),] 

#remove variable column
long_GO <- long_GO[, c(1, 3)]

#sort by transcript/gene
gene.go <- long_GO[order(long_GO$V1), ]

#create list with element for each gene, containing vectors with all terms for each gene
gene2GO <- tapply(gene.go$`IDs`, gene.go$V1, function(x)x)

head(gene2GO) #go IDs as strings

background <- colnames(mat)

#col= 'turquoise'
col='blue'


geneList <- colnames(mat[,ModuleColors==col]) ##this must get the ones
NOTgene_list <- setdiff(background, geneList) ##this must get the zeros

temp <- c(rep(1, each = length(geneList)), rep(0, each = length(NOTgene_list)))
names(temp) <- c(geneList, NOTgene_list)

                  
#Create topGOdata object; rerun this code if you want to explore different ontologies
GOdata <- new("topGOdata",
                   ontology = "BP",  #ontology criteria = biological process
                   allGenes = temp, #gene universe because all genes are present in the list; here only the ones with read abundance threshold are included because DESeq only worked with them
                   geneSelectionFun = function(x)(x == 1), ##function allows to use the genes which are de for treatment
                   annot = annFUN.gene2GO, gene2GO = gene2GO) #gene ID-to-GO terms

#run enrichment test: here Fishers Exact Test and explore how the algorithms perform differently to decorrelate the GO graph structure
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
resultFisher.weight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

#summarize in table
GOtable <- GenTable(GOdata, p.value.elim = resultFisher.elim, p.value.weight01 = resultFisher.weight01, topNodes = length(resultFisher.weight01@score), numChar = 10000000)

write.csv(GOtable, file = paste(col, '_WGCNA_module_BP_enrichment_fisher.csv', sep = ''))
#write.csv(GOtable, file = paste(col, '_WGCNA_module_MF_enrichment_fisher.csv', sep = ''))

gg1
dev.off()
