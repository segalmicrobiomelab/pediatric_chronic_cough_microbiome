##Load the files needed
file = "otu_table.Merged.BIOLandBKG.f.s.biom"

#includes RI and RI cut-offs of 4 and 6 
map = "Map.Pediatric.Merged.final.B1.txt"

#install the required packages
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("d3heatmap")
install.packages("vegan")
install.packages("Heatplus")
install.packages("igraph")
source("https://bioconductor.org/biocLite.R")
biocLite("Heatplus")
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
biocLite("stringi", type = "source")
biocLite("phyloseq")
install.packages('ggplot2')
install.packages("ape")
install.packages("ade4")
install.packages("phyloseq")
install.packages("ade4")
install.packages("stringi")


#Load Phyloseq
library(phyloseq)
#Load GG Plot 2
library(ggplot2)

load(file="180103.Pediatric.RData")
save.image(file="180103.Pediatric.RData")

#Load other libraries: 
library(ade4)
library("RColorBrewer")
library('ape')
library("plyr")
library("gplots")
library("d3heatmap")
library("vegan")
library("vegan3d")
library("Heatplus")
library("igraph")
library("reshape2")
library("MASS")

theme_set(theme_bw())




# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))

lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)


#Give a colnames to separate different taxonomic levels

colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Load the tree file (use the unannotated.tree)
treefile = "97_otus_unannotated.tree"
tree.obj = import_qiime(treefilename = treefile) 


# Now merge the three separate phyloseq objects into a single object
otu.table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)

rownames(sample_data(otu.table))
colnames(sample_data(otu.table))


# Remove taxa with 0 abundance
otu.table = subset_taxa(otu.table, rowSums(otu_table(otu.table)) != 0)


##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)

colnames(sample_data(otu.relative.table))
rownames(sample_data(otu.relative.table))


# Create phyllum and order tables (do it after normalization and out of the relative table)
Phylum.rel.table = tax_glom(otu.relative.table, taxrank = "Phylum")
Class.rel.table = tax_glom(otu.relative.table, taxrank = "Class")
Order.rel.table = tax_glom(otu.relative.table, taxrank = "Order")
Family.rel.table = tax_glom(otu.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(otu.relative.table, taxrank = "Genus")
OTU.rel.table = tax_glom(otu.relative.table, taxrank = "OTU") 

colnames(sample_data(OTU.rel.table))



########################################
#Figure 1

#######HEAT MAPS
#install necessary packages - one time only 
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("d3heatmap")
install.packages("vegan")

#load libraries 
library("gplots")
library("RColorBrewer")
library("d3heatmap")
library("vegan")
library("devtools")


##Set of pruning we used:
##Here we can select for genera present in >5% relative abundance in 0.5% of the samples (this approach brings in > 70% of the data in almost all samples)
Genus.Rel.wh1 = genefilter_sample(Genus.rel.table, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(Genus.rel.table))
Genus.Rel.table1B = prune_taxa(Genus.Rel.wh1, Genus.rel.table)
colnames(sample_data(Genus.Rel.table1B)) ##will show metadata columns
rownames(sample_data(Genus.Rel.table1B)) ##will show metadata columns

pdf("Genus Pruned Bar Plot.pdf", width = 28, height = 8)
plot_bar(Genus.Rel.table1B, fill="Genus")
dev.off()

#set data tables  
data <- otu_table(OTU.rel.table) #all data 
GenusData <-otu_table(Genus.Rel.table1B) #pruned to selected Genuses based on abundance

#change colors to 0 = white, 1 = blue 
mypalette <- colorRampPalette(c('#ffffff','#4169E1','#0000CD'))


#create vector to lable by Topograph Code
TopoVector = sample_data(OTU.rel.table)$Topograph_code
#duplicate to create a color vector and replace value w/ color 
#Colorvector can only replace numbers! 
Colorvector <-TopoVector
Colorvector <- replace(Colorvector, which (Colorvector == "1"), c("green")) #BKG
Colorvector <- replace(Colorvector, which (Colorvector == "2"), "red") #Upper A
Colorvector <- replace(Colorvector, which (Colorvector == "3"), "blue") #lower A
Colorvector <- replace(Colorvector, which (Colorvector == "4"), "cadetblue1") #Middle A
Colorvector <- replace(Colorvector, which (Colorvector == "5"), "orange") #upper GI
Colorvector <- replace(Colorvector, which (Colorvector == "6"), "purple") #Middle GI
Colorvector <- replace(Colorvector, which (Colorvector == "7"), "darkmagenta") #lower GI


#Here we are able to change the names for genuses that are labelled as "g__" --> Come back to this
x10 = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.table1B))), ntaxa(Genus.Rel.table1B)), Genus.Rel.table1B)
tax_table(Genus.Rel.table1B)
# Add a new rank, Strain, with the Genus ids
tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
# Define the ranks you want to include
myranks = c("Order", "Family", "Genus")
mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")
# Add concatenated labels as a new rank after strain
tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)
# Check this out on a tree
plot_tree(x10, label.tips="catglab", color="topograph", nodelabf=nodeplotboot(),ladderize="left", plot.margin=2.15, size="abundance")


##Cluster Heatmap

#cluster Genuses(row) by Bray 
GenusData.Bray.dist <-vegdist(GenusData, method = "bray")
Genus.Bray.clus <-hclust(GenusData.Bray.dist, "aver")

#cluster genuses by unifrac 
GenusData.wUnifrac.dist =distance(Genus.Rel.table1B, method = "wUnifrac")
Genus.wUnifrac.clus <-hclust(GenusData.wUnifrac.dist, "aver")


#cluster samples(Col) by Bray
Samples.Bray.dist = distance(GenusData, method="bray")
Samples.cluster.Bray = hclust(Samples.Bray.dist, "aver")

#cluster samples(Col) by Unifrac
Samples.wUnifrac.dist = distance(Genus.Rel.table1B, method="wunifrac", type="samples")
Samples.cluster.wUnifrac = hclust(Samples.wUnifrac.dist, "aver")



## TO PLOT

#Bray by unifrac
pdf("Bray_Unifrac_Heatmap_AllSamples.pdf", height = 10, width = 20)
heatmap.2(GenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(Genus.Rel.table1B)$Paper_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance",
)
dev.off()


#Subset by type from Genus.Rel.table1B
BKG.Rel.Table = subset_samples(Genus.Rel.table1B,Topograph_code==1)
UA.Rel.Table = subset_samples(Genus.Rel.table1B, Topograph_code==2)
LA.Rel.Table = subset_samples(Genus.Rel.table1B, Topograph_code==3)
MA.Rel.Table = subset_samples(Genus.Rel.table1B, Topograph_code==4)
UGI.Rel.Table = subset_samples(Genus.Rel.table1B, Topograph_code==5)
MGI.Rel.Table = subset_samples(Genus.Rel.table1B, Topograph_code==6)
LGI.Rel.Table = subset_samples(Genus.Rel.table1B, Topograph_code==7)

#cluster samples types by Unifrac 
BKG.Samples.wUnifrac.dist = distance(BKG.Rel.Table, method="wunifrac", type="samples")
BKG.Samples.cluster.wUnifrac = hclust(BKG.Samples.wUnifrac.dist, "aver")

UA.Samples.wUnifrac.dist = distance(UA.Rel.Table, method="wunifrac", type="samples")
UA.Samples.cluster.wUnifrac = hclust(UA.Samples.wUnifrac.dist, "aver")

LA.Samples.wUnifrac.dist = distance(LA.Rel.Table, method="wunifrac", type="samples")
LA.Samples.cluster.wUnifrac = hclust(LA.Samples.wUnifrac.dist, "aver")

MA.Samples.wUnifrac.dist = distance(MA.Rel.Table, method="wunifrac", type="samples")
MA.Samples.cluster.wUnifrac = hclust(MA.Samples.wUnifrac.dist, "aver")

UGI.Samples.wUnifrac.dist = distance(UGI.Rel.Table, method="wunifrac", type="samples")
UGI.Samples.cluster.wUnifrac = hclust(UGI.Samples.wUnifrac.dist, "aver")

MGI.Samples.wUnifrac.dist = distance(MGI.Rel.Table, method="wunifrac", type="samples")
MGI.Samples.cluster.wUnifrac = hclust(MGI.Samples.wUnifrac.dist, "aver")

LGI.Samples.wUnifrac.dist = distance(LGI.Rel.Table, method="wunifrac", type="samples")
LGI.Samples.cluster.wUnifrac = hclust(LGI.Samples.wUnifrac.dist, "aver")

sample_data(Genus.Rel.table1B)$Topograph_code
sample_data(Genus.Rel.table1B)$SampleType

#Genus OTU Data must also be separated 
BKGGenusData <-otu_table(BKG.Rel.Table)
UAGenusData <-otu_table(UA.Rel.Table)
LAGenusData <-otu_table(LA.Rel.Table)
MAGenusData <-otu_table(MA.Rel.Table)
UGIGenusData <-otu_table(UGI.Rel.Table)
MGIGenusData <-otu_table(MGI.Rel.Table)
LGIGenusData <-otu_table(LGI.Rel.Table)


#create heatmaps of each subset w/ same bray distance calculated taxa 
#BKG
pdf("Bray_Unifrac_Heatmap_BKG.a.pdf", height = 10, width = 20)
heatmap.2(BKGGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(BKG.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(BKG.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Background",
)
dev.off()

#UA
pdf("Bray_Unifrac_Heatmap_UA.a3.pdf", height = 10, width = 10)
heatmap.2(UAGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(UA.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(UA.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Upper Airway",
)
dev.off()

#LA
pdf("Bray_Unifrac_Heatmap_LA.a3.pdf", height = 10, width = 10)
heatmap.2(LAGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(LA.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(LA.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Lower Airway",
)
dev.off()

#MA 
pdf("Bray_Unifrac_Heatmap_MA.a3.pdf", height = 10, width = 10)
heatmap.2(MAGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(MA.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(MA.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Middle Airway",
)
dev.off()

#UGI
pdf("Bray_Unifrac_Heatmap_UGI.a3.pdf", height = 10, width = 10)
heatmap.2(UGIGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(UGI.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(UGI.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Upper GI",
)
dev.off()

#MGI
pdf("Bray_Unifrac_Heatmap_MGI.a3.pdf", height = 10, width = 10)
heatmap.2(MGIGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(MGI.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(MGI.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Middle GI",
)
dev.off()

#LGI
pdf("Bray_Unifrac_Heatmap_LGI.a3.pdf", height = 10, width = 10)
heatmap.2(LGIGenusData, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(LGI.Samples.cluster.wUnifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(LGI.Rel.Table)$Subject_ID,
	cexCol = .35,
	col = mypalette(17),
 	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks =c(seq(0,.1,length=10),seq(.11,0.3,length=4),seq(0.31,.7,length=4)),
	#ColSideColors=Colorvector,
	main = "Heatmap of Unifrac Distance - Lower GI",
)
dev.off()



################################################
#Figure 2 Panel A
## Alpha Diversity =  
Shannon_diversity = diversity(otu_table(otu.relative.table), index = "shannon", MARGIN = 2, base = exp(1))

##First check figure before relabeling
boxplot(Shannon_diversity ~ sample_data(otu.relative.table)$topograph)

#Reorder factors
sample_data(otu.relative.table)$topograph<-factor(sample_data(otu.relative.table)$topograph, levels = c("BKG", "U.A", "M.A", "L.A", "U.GI", "M.GI", "L.GI"))
boxplot(Shannon_diversity ~ sample_data(otu.relative.table)$topograph)

## TO PLOT IT
pdf(file="Shannon_otu_diversity_labels.pdf", width=12, height=4)
boxplot(Shannon_diversity ~ sample_data(otu.relative.table)$topograph)
dev.off()


#Significance testing
kruskal.test(Shannon_diversity ~ sample_data(otu.relative.table)$topograph)
#Result
	Kruskal-Wallis rank sum test

data:  Shannon_diversity by sample_data(otu.relative.table)$topograph
Kruskal-Wallis chi-squared = 68.445, df = 6, p-value = 8.517e-13

#Table used to create figure in Prism
write.table(Shannon_diversity, file="Shannon.diversity.txt", sep="\t")





################################################
# Figure 2 Panel B
####To calculate W Unifrac Distance
library(ade4)
#Calculate Distance
wUniF.dist = UniFrac(otu.relative.table, weighted=TRUE)
#estimate number of axes
wUniF.pco = dudi.pco(cailliez(wUniF.dist))

##Plot PCoA
pdf(file="Pediatric.UniF.a.pdf", width=16, height=16)
s.class(wUniF.pco $li, sample_data(otu.relative.table)$topograph, label=c("Background", "BAL", "Duodenum", "Trachea", "Gastric", "Supraglottic", "Esophagus"), col=c("green", "blue", "darkmagenta ", "cadetblue1 ", "purple", "red", "orange"))
dev.off()

pdf(file="Pediatric.UniF.b.pdf", width=16, height=16)
s.class(wUniF.pco $li, sample_data(otu.relative.table)$topograph, label=c("Background", "BAL", "Duodenum", "Trachea", "Gastric", "Tech.Control", "Supraglottic", "Esophagus"), col=c("green", "blue", "darkmagenta ", "cadetblue1 ", "purple", "black", "red", "orange"))
dev.off()


#calculate difference
adonis(wUniF.dist ~ topograph, data=data.frame(sample_data(otu.relative.table)))
#Result:
Call:
adonis(formula = wUniF.dist ~ topograph, data = data.frame(sample_data(otu.relative.table))) 

Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
topograph   6    3.6831 0.61385  9.0079 0.13477  0.001 ***
Residuals 347   23.6467 0.06815         0.86523           
Total     353   27.3298                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


ordinate.wunifrac <- ordinate(otu.relative.table, method="PCoA", distance="wunifrac")

plot_ordination(otu.relative.table, ordinate.wunifrac, color="topograph", shape="topograph")

########################################
#Figure 4 Panel A
## Alpha Diversity  
Shannon_diversity_BAL_Classification = diversity(otu_table(BAL.otu.relative.table), index = "shannon", MARGIN = 2, base = exp(1))
Sipmson_diversity_BAL_Classification = diversity(otu_table(BAL.otu.relative.table), index = "simpson", MARGIN = 2, base = exp(1))


##First check figure before relabeling
boxplot(Shannon_diversity_BAL_Classification ~ sample_data(BAL.otu.relative.table)$Classification)

#Reorder factors
sample_data(BAL.otu.relative.table)$Classification <-factor(sample_data(BAL.otu.relative.table)$Classification, levels = c("Asthma", "Bacterial.Bronchitis", "Neuro.Impaired.PO", "Neuro.Impaired.Enteral"))
boxplot(Shannon_diversity_BAL_Classification ~ sample_data(BAL.otu.relative.table)$Classification)



##Plot Alpha Diversity Boxplot
pdf(file="Shannon_Classifcation_BAL_otu_diversity_labels.pdf", width=6, height=12)
boxplot(Shannon_diversity_BAL_Classification ~ sample_data(BAL.otu.relative.table)$Classification, col = c("darkgreen", "blue", "red", "darkred"), names  = c("Asthma", "Bacterial Bronchitis", "Neuro Impaired PO Feed", "Neuro Impaired Enteral Feed"))
dev.off()

#Significance testing
kruskal.test(Shannon_diversity_BAL_Classification ~ sample_data(BAL.otu.relative.table)$Classification)
##Result
	Kruskal-Wallis rank sum test

data:  Shannon_diversity_BAL_Classification by sample_data(BAL.otu.relative.table)$Classification
Kruskal-Wallis chi-squared = 8.8125, df = 3, p-value = 0.03189

#Output table used to make boxplot in Prism
write.table(Shannon_diversity_BAL_Classification, file="Shannon.diversity.BAL.Classification.txt", sep="\t")




########################################
#Figure 4 Panel B
#Beta Diversity
BAL.otu.relative.table = subset_samples(OTU.rel.table, Topograph_code==3) 
colnames(sample_data(BAL.otu.relative.table))
#Calculate Distance
BAL.wUniF.dist = UniFrac(BAL.otu.relative.table, weighted=TRUE)
#estimate number of axes
BAL.wUniF.pco = dudi.pco(cailliez(BAL.wUniF.dist))

colnames(sample_data(BAL.otu.relative.table))
rownames(sample_data(BAL.otu.relative.table))

#Plot first to see colors
s.class(BAL.wUniF.pco$li, sample_data(BAL.otu.relative.table)$Classification, col=c("darkgreen", "blue", "darkred", "red"))

#Relabel PCoA
s.class(BAL.wUniF.pco$li, sample_data(BAL.otu.relative.table)$Classification, label=c("Asthma", "Bacterial Bronchitis", "Neuro Impaired PO Fed", "Neuro Impaired Enteral Fed"), col=c("darkgreen", "blue", "darkred", "red"))


pdf(file="Pediatric.BAL.Classification.UniF.withRI.pdf", width=12, height=12)
s.class(BAL.wUniF.pco$li, sample_data(BAL.otu.relative.table)$Classification, label=c("Asthma", "Bacterial Bronchitis", "Neuro Impaired PO Fed", "Neuro Impaired Enteral Fed"), col=c("darkgreen", "blue", "darkred", "red"))
dev.off()


adonis(BAL.wUniF.dist ~ Classification, data=data.frame(sample_data(BAL.otu.relative.table)))
#Result
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
Classification  3   0.45162 0.150539  2.4677 0.18788  0.008 **
Residuals      32   1.95211 0.061004         0.81212          
Total          35   2.40373                  1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
 
################### Figure 4 Panel B separated by RI > 4 data 

#Plot by RI > 4
colnames(sample_data(BAL.otu.relative.table))
#remove those without an RI score
#BAL.otu.relative.table.RI = subset_samples(BAL.otu.relative.table, RI_impedence_timereflux !="n.a")

s.class(BAL.wUniF.pco$li, sample_data(BAL.otu.relative.table)$RI_4,  col=c("darkgreen", "blue", "red"))
label=c("RI < 4", "RI > 4"),
adonis(BAL.wUniF.dist ~ Classification, data=data.frame(sample_data(BAL.otu.relative.table)))
#Result
Permutation: free
Number of permutations: 999

Terms added sequentially (first to last)

               Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
Classification  3   0.45941 0.153137  2.4469 0.18659  0.005 **
Residuals      32   2.00268 0.062584         0.81341          
Total          35   2.46209                  1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#Plot by RI>6

s.class(BAL.wUniF.pco$li, sample_data(BAL.otu.relative.table)$RI_6,  col=c("darkgreen", "blue", "red"))
label=c("RI < 4", "RI > 4"),



########################################
#Figure 6
###Biplots for Classification by wUniFrac
##Here we can select for genera present in >5% relative abundance in 0.5% of the samples (this approach brings in > 75% of the data in almost all samples)
LA.Genus.Rel.table = subset_samples(Genus.rel.table, Topograph_code==3)
rownames(sample_data(LA.Genus.Rel.table))
LA.Genus.Rel.wh1 = genefilter_sample(LA.Genus.Rel.table, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(LA.Genus.Rel.table))
LA.Genus.Rel.table1B = prune_taxa(LA.Genus.Rel.wh1, LA.Genus.Rel.table)
colnames(sample_data(LA.Genus.Rel.table1B))
plot_bar(LA.Genus.Rel.table1B, fill="Genus")

LA.Genus.ords <- ordinate(LA.Genus.Rel.table1B, "NMDS", "wunifrac")
plot_ordination(LA.Genus.Rel.table1B, LA.Genus.ords, type = "biplot", color = "Genus", title = "Lower Airway Biplot", shape = "Classification", label="Subject_ID") + geom_point(size = 2)

colnames(sample_data(LA.Genus.ords))


#Use this file with taxa in different colors so you can edit the figure in illustrator accurately
pdf("BAL_Biplot.Classification.wUnifrac.pdf", width = 12, height = 12)
plot_ordination(LA.Genus.Rel.table1B, LA.Genus.ords, type = "biplot", color = "Genus", title = "Lower Airway Biplot", shape = "Classification", label= "Genus") + geom_point(size = 2) + scale_shape_manual(values =c(16, 17, 18, 19, 15))
dev.off()

#For labeling unknown general
pdf("BAL_Biplot.Classification.family.wUnifrac.pdf", width = 12, height = 12)
plot_ordination(LA.Genus.Rel.table1B, LA.Genus.ords, type = "biplot", color = "Genus", title = "Lower Airway Biplot", shape = "Classification", label= "Family") + geom_point(size = 2) + scale_shape_manual(values =c(16, 17, 18, 19, 15))
dev.off()

#Final figure, needs to be edited in Illustrator
pdf("BAL_Biplot.Classification.a.wUnifrac.pdf", width = 12, height = 12)
plot_ordination(LA.Genus.Rel.table1B, LA.Genus.ords, type = "biplot", color = "Classification", title = "Lower Airway Biplot", shape = "Classification", label= "Genus") + geom_point(size = 2) + scale_shape_manual(values =c(16, 17, 18, 19, 15)) + scale_color_manual(values = c("#F45E5A", "#006B01", "#0201FE", "#800002", "#FF0000"))
dev.off()



########################################
#Figure 6 - Vectorgram Overlay
#Note: This code was used to create the vectorgram by wUniFrac, but used pruned data. The code used to create the stats in Table 2 can be found further below.

#subset only samples for BAL 
BAL.Genus.Rel.table = subset_samples(Genus.rel.table, Topograph_code==3)


##Here we can select for genera present in >5% relative abundance in 0.5% of the samples (this approach brings in > 75% of the data in almost all samples)
rownames(sample_data(BAL.Genus.Rel.table))
BAL.Genus.Rel.wh1 = genefilter_sample(BAL.Genus.Rel.table, filterfun_sample(function(x) x > 0.05), A = 0.005 * nsamples(BAL.Genus.Rel.table))
BAL.Genus.Rel.table1B = prune_taxa(BAL.Genus.Rel.wh1, BAL.Genus.Rel.table)
colnames(sample_data(BAL.Genus.Rel.table1B))
plot_bar(BAL.Genus.Rel.table1B, fill="Genus")

#subset only those w/ RI scores
BAL.Genus.Rel.RI.table1B = subset_samples(BAL.Genus.Rel.table1B, RI_impedence_timereflux !="n.a")


#Unifrac
BAL.Unifrac.filt.pruned.RI.dist = distance(BAL.Genus.Rel.RI.table1B, "wunifrac")


######Unifrac 
BAL.Unifrac.pruned.RI.mds <- metaMDS(BAL.Unifrac.filt.pruned.RI.dist, trace = FALSE)
BAL.Unifrac.pruned.RI.mds
Call:
metaMDS(comm = BAL.Unifrac.filt.pruned.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     BAL.Unifrac.filt.pruned.RI.dist 
Distance: user supplied 

Dimensions: 2 
Stress:     0.1468912 
Stress type 1, weak ties
No convergent solutions - best solution after 20 tries
Scaling: centring, PC rotation 
Species: scores missing


ordiplot(BAL.Unifrac.pruned.RI.mds, type = "t")
ordiplot(BAL.Unifrac.pruned.RI.mds, type = "p")

BAL.Classification.RI.map = "BAL.mv.metadata.continuous.RI.txt"
BAL.RI.pruned.table=sample_data(read.table(BAL.Classification.RI.map, header=T, sep="\t", row.names=1))
BAL.RI.pruned.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!

BAL.Unifrac.RI.pruned.ef <- envfit(BAL.Unifrac.pruned.RI.mds, BAL.RI.pruned.table, permu = 999)
BAL.Unifrac.RI.pruned.ef

plot(BAL.Unifrac.pruned.RI.mds, type = "t", display = "sites")
plot(BAL.Unifrac.RI.pruned.ef, p.max=0.999)
***VECTORS

                           NMDS1    NMDS2     r2 Pr(>r)   
BSS_revised              0.91231  0.40949 0.3103  0.010 **
Macrophages             -0.94612 -0.32381 0.3227  0.013 * 
Lymphocytes             -0.87188 -0.48971 0.0406  0.591   
Neutrophils              0.86181  0.50722 0.2238  0.048 * 
Eosinophils              0.29630  0.95509 0.0313  0.690   
Total_cells              0.44131  0.89735 0.1007  0.261   
RI_impedence_timereflux  0.99945  0.03327 0.0100  0.875   
age_today_dob           -0.67647  0.73647 0.0427  0.620   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

pdf(file = "BAL.Unifrac.Vectorgram.pruned.RI.pdf", width = 12, height = 12)
plot(BAL.Unifrac.pruned.RI.mds, type = "t", display = "sites")
plot(BAL.Unifrac.RI.pruned.ef, p.max=0.999)
dev.off()

pdf(file= "BAL.Unifrac.VectorgramONLY.pruned.RI.pdf", width = 12, height = 12)
plot(BAL.Unifrac.pruned.RI.mds, type = "n")
plot(BAL.Unifrac.RI.pruned.ef, p.max=0.999)
dev.off()


########################################
#Supplementary Table 1 (Vectorgram statistical testing)
#BAL
#NOT PRUNED BAL vectorgram for RI score data
#subset only samples with data for pH composite score

BAL.otu.relative.table.RI = subset_samples(BAL.otu.relative.table, RI_impedence_timereflux !="n.a")

rownames(sample_data(BAL.otu.relative.table.RI))

#Bray distance
BAL.Bray.RI.dist = distance(BAL.otu.relative.table.RI, "bray")

BAL.Bray.RI.mds <- metaMDS(BAL.Bray.RI.dist, trace = FALSE)
BAL.Bray.RI.mds

#Call:
metaMDS(comm = BAL.Bray.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     BAL.Bray.RI.dist 
Distance: bray 

Dimensions: 2 
Stress:     0.1336178 
Stress type 1, weak ties
Two convergent solutions found after 20 tries
Scaling: centring, PC rotation 
Species: scores missing

ordiplot(BAL.Bray.RI.mds, type = "t")
ordiplot(BAL.Bray.RI.mds, type = "p")


BAL.Classification.RI.map = "BAL.mv.metadata.continuous.RI.txt"
BAL.RI.table=sample_data(read.table(BAL.Classification.RI.map, header=T, sep="\t", row.names=1))
BAL.RI.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!


BAL.Bray.RI.ef <- envfit(BAL.Bray.RI.mds, BAL.RI.table, permu = 999)
BAL.Bray.RI.ef

#Result
***VECTORS

                           NMDS1    NMDS2     r2 Pr(>r)  
BSS_revised             -0.58974 -0.80759 0.3079  0.017 *
Macrophages              0.87117 -0.49097 0.0404  0.608  
Lymphocytes              0.37349  0.92764 0.1844  0.089 .
Neutrophils             -0.45120 -0.89242 0.0747  0.421  
Eosinophils              0.08626 -0.99627 0.0943  0.308  
Total_cells             -0.05897 -0.99826 0.1466  0.136  
RI_impedence_timereflux -0.44695  0.89456 0.0862  0.340  
age_today_dob            0.75799 -0.65227 0.0157  0.811  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

plot(BAL.Bray.RI.mds, type = "t", display = "sites")
plot(BAL.Bray.RI.ef, p.max=0.999)


pdf(file = "BAL.Bray.Vectorgram.RI.pdf", width = 12, height = 12)
plot(BAL.Bray.pH.mds, type = "t", display = "sites")
plot(BAL.Bray.pH.ef, p.max=0.999)
dev.off()


######Supraglottic
#NOT PRUNED Sup vectorgram for RI score data
#subset to Sup only
Sup.otu.relative.table = subset_samples(OTU.rel.table, SampleType=="Sup") 

#subset only samples with data for RI composite score

Sup.otu.relative.table.RI = subset_samples(Sup.otu.relative.table, RI_impedence_timereflux !="n.a")

rownames(sample_data(Sup.otu.relative.table.RI))

#Bray distance
Sup.Bray.RI.dist = distance(Sup.otu.relative.table.RI, "bray")


Sup.Bray.RI.mds <- metaMDS(Sup.Bray.RI.dist, trace = FALSE)
Sup.Bray.RI.mds

#Call:
metaMDS(comm = Sup.Bray.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     Sup.Bray.RI.dist 
Distance: bray 

Dimensions: 2 
Stress:     0.1391274 
Stress type 1, weak ties
Two convergent solutions found after 20 tries
Scaling: centring, PC rotation 
Species: scores missing



ordiplot(Sup.Bray.RI.mds, type = "t")
ordiplot(Sup.Bray.RI.mds, type = "p")

Sup.Classification.RI.map = "Sup.mv.metadata.continuous.RI.txt"
Sup.RI.table=sample_data(read.table(Sup.Classification.RI.map, header=T, sep="\t", row.names=1))
Sup.RI.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!

Sup.Bray.RI.ef <- envfit(Sup.Bray.RI.mds, Sup.RI.table, permu = 999)
Sup.Bray.RI.ef


#Result
***VECTORS

                           NMDS1    NMDS2     r2 Pr(>r)
BSS_revised             -0.86461 -0.50245 0.0271  0.712
Macrophages              0.27448  0.96159 0.0015  0.978
Lymphocytes              0.48883  0.87238 0.0368  0.583
Neutrophils             -0.65245 -0.75783 0.0165  0.817
Eosinophils              0.30062  0.95375 0.0660  0.417
Total_cells              0.08137 -0.99668 0.0074  0.830
age_today_dob            0.25770 -0.96622 0.0415  0.609
RI_impedence_timereflux  0.11626 -0.99322 0.0488  0.519

---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

plot(Sup.Bray.RI.mds, type = "t", display = "sites")
plot(Sup.Bray.RI.ef, p.max=0.999)


pdf(file = "Sup.Bray.Vectorgram.RI.pdf", width = 12, height = 12)
plot(Sup.Bray.RI.mds, type = "t", display = "sites")
plot(Sup.Bray.RI.ef, p.max=0.999)
dev.off()




#Trachea
##NOT PRUNED Trach vectorgram for RI score data
#subset to Trach only
Trach.otu.relative.table = subset_samples(OTU.rel.table, Topograph_code==4) 

colnames(sample_data(OTU.rel.table))

#subset only samples with data for RI 
Trach.otu.relative.table.RI = subset_samples(Trach.otu.relative.table, RI_impedence_timereflux !="n.a")

rownames(sample_data(Trach.otu.relative.table.RI))

#Bray distance
Trach.Bray.RI.dist = distance(Trach.otu.relative.table.RI, "bray")


Trach.Bray.RI.mds <- metaMDS(Trach.Bray.RI.dist, trace = FALSE)
Trach.Bray.RI.mds

#Call:
metaMDS(comm = Trach.Bray.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     Trach.Bray.RI.dist 
Distance: bray 

Dimensions: 2 
Stress:     0.1888605 
Stress type 1, weak ties
No convergent solutions - best solution after 20 tries
Scaling: centring, PC rotation 
Species: scores missing

ordiplot(Trach.Bray.RI.mds, type = "t")
ordiplot(Trach.Bray.RI.mds, type = "p")


Trach.Classification.RI.map = "Trach.mv.metadata.continuous.RI.txt"
Trach.RI.table=sample_data(read.table(Trach.Classification.RI.map, header=T, sep="\t", row.names=1))
Trach.RI.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!


Trach.Bray.RI.ef <- envfit(Trach.Bray.RI.mds, Trach.RI.table, permu = 999)
Trach.Bray.RI.ef

#Results
***VECTORS

                      NMDS1    NMDS2     r2 Pr(>r)   
BSS_revised             -0.94736  0.32017 0.3843  0.001 ***
Macrophages              0.99690 -0.07868 0.3285  0.010 ** 
Lymphocytes              0.99819 -0.06008 0.1073  0.227    
Neutrophils             -0.99810  0.06164 0.4201  0.002 ** 
Eosinophils             -0.58503  0.81101 0.0593  0.492    
Total_cells             -0.44279  0.89663 0.1333  0.126    
age_today_dob           -0.39244  0.91978 0.0051  0.939    
RI_impedence_timereflux -0.91027  0.41402 0.0113  0.872    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

plot(Trach.Bray.RI.mds, type = "t", display = "sites")
plot(Trach.Bray.RI.ef, p.max=0.999)


pdf(file = "Trach.Bray.Vectorgram.RI.pdf", width = 12, height = 12)
plot(Trach.Bray.RI.mds, type = "t", display = "sites")
plot(Trach.Bray.RI.ef, p.max=0.999)
dev.off()


#Esophagus
#NOT PRUNED Esoph vectorgram for RI composite score data
#subset to Esoph only
Esoph.otu.relative.table = subset_samples(OTU.rel.table, Topograph_code==5) 
sample_data(OTU.rel.table)$SampleType

#subset only samples with data for pH composite score
Esoph.otu.relative.table.RI = subset_samples(Esoph.otu.relative.table,  RI_impedence_timereflux!="n.a")

rownames(sample_data(Esoph.otu.relative.table.RI))

#Bray distance
Esoph.Bray.RI.dist = distance(Esoph.otu.relative.table.RI, "bray")

Esoph.Bray.RI.mds <- metaMDS(Esoph.Bray.RI.dist, trace = FALSE)
Esoph.Bray.RI.mds

#Call:
metaMDS(comm = Esoph.Bray.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     Esoph.Bray.RI.dist 
Distance: bray 

Dimensions: 2 
Stress:     0.1561063 
Stress type 1, weak ties
Two convergent solutions found after 20 tries
Scaling: centring, PC rotation 
Species: scores missing

ordiplot(Esoph.Bray.RI.mds, type = "t")
ordiplot(Esoph.Bray.RI.mds, type = "p")


Esoph.Classification.RI.map = "Esoph.mv.metadata.continuous.RI.txt"
Esoph.RI.table=sample_data(read.table(Esoph.Classification.RI.map, header=T, sep="\t", row.names=1))
Esoph.RI.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!

Esoph.Bray.RI.ef <- envfit(Esoph.Bray.RI.mds, Esoph.RI.table, permu = 999)
Esoph.Bray.RI.ef


#Result
                           NMDS1    NMDS2     r2 Pr(>r)
BSS_revised             -0.95337  0.30182 0.0793  0.363
Macrophages              0.96906 -0.24681 0.0424  0.592
Lymphocytes              0.23955  0.97088 0.1384  0.166
Neutrophils             -0.57519 -0.81802 0.0700  0.445
Eosinophils             -0.99780  0.06626 0.0381  0.627
Total_cells             -0.59353  0.80481 0.0125  0.759
age_today_dob            0.99902  0.04427 0.0991  0.273
RI_impedence_timereflux  0.97093 -0.23937 0.0555  0.452
 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

plot(Esoph.Bray.RI.mds, type = "t", display = "sites")
plot(Esoph.Bray.RI.ef, p.max=0.999)


pdf(file = "Esoph.Bray.Vectorgram.RI.pdf", width = 12, height = 12)
plot(Esoph.Bray.RI.mds, type = "t", display = "sites")
plot(Esoph.Bray.RI.ef, p.max=0.999)
dev.off()




#Gastric
#NOT PRUNED Gast vectorgram for RI composite score data

Gast.otu.relative.table = subset_samples(OTU.rel.table, Topograph_code==6)
#subset only samples with data for RI composite score
Gast.otu.relative.table.RI = subset_samples(Gast.otu.relative.table, RI_impedence_timereflux !="n.a")

rownames(sample_data(Gast.otu.relative.table.RI))

#Bray distance
Gast.Bray.RI.dist = distance(Gast.otu.relative.table.RI, "bray")


Gast.Bray.RI.mds <- metaMDS(Gast.Bray.RI.dist, trace = FALSE)
Gast.Bray.RI.mds

#Call:
metaMDS(comm = Gast.Bray.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     Gast.Bray.RI.dist 
Distance: bray 

Dimensions: 2 
Stress:     0.1806728 
Stress type 1, weak ties
No convergent solutions - best solution after 20 tries
Scaling: centring, PC rotation 
Species: scores missing



ordiplot(Gast.Bray.RI.mds, type = "t")
ordiplot(Gast.Bray.RI.mds, type = "p")


Gast.Classification.RI.map = "Gast.mv.metadata.continuous.RI.txt"
Gast.RI.table=sample_data(read.table(Gast.Classification.RI.map, header=T, sep="\t", row.names=1))
Gast.RI.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!

Gast.Bray.RI.ef <- envfit(Gast.Bray.RI.mds, Gast.RI.table, permu = 999)
Gast.Bray.RI.ef


#Result 
***VECTORS

                           NMDS1    NMDS2     r2 Pr(>r)
BSS_revised              0.25216 -0.96769 0.0458  0.587
Macrophages             -0.99936 -0.03567 0.0100  0.872
Lymphocytes             -0.63772  0.77027 0.1470  0.147
Neutrophils              0.81045 -0.58580 0.0769  0.359
Eosinophils             -0.87370  0.48647 0.0514  0.527
Total_cells              0.52118 -0.85345 0.0090  0.826
age_today_dob           -0.89228 -0.45149 0.1445  0.166
RI_impedence_timereflux  0.21074 -0.97754 0.0114  0.843  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

plot(Gast.Bray.RI.mds, type = "t", display = "sites")
plot(Gast.Bray.RI.ef, p.max=0.999)


pdf(file = "Gast.Bray.Vectorgram.RI.pdf", width = 12, height = 12)
plot(Gast.Bray.RI.mds, type = "t", display = "sites")
plot(Gast.Bray.RI.ef, p.max=0.999)
dev.off()


#Duodenum
#NOT PRUNED Duod vectorgram for RI composite score data

Duod.otu.relative.table = subset_samples(OTU.rel.table, Topograph_code==7)
#subset only samples with data for RI composite score
Duod.otu.relative.table.RI = subset_samples(Duod.otu.relative.table, RI_impedence_timereflux !="n.a")

rownames(sample_data(Duod.otu.relative.table.RI))

#Bray distance
Duod.Bray.RI.dist = distance(Duod.otu.relative.table.RI, "bray")


Duod.Bray.RI.mds <- metaMDS(Duod.Bray.RI.dist, trace = FALSE)
Duod.Bray.RI.mds

#Call:
metaMDS(comm = Duod.Bray.RI.dist, trace = FALSE) 

global Multidimensional Scaling using monoMDS

Data:     Duod.Bray.RI.dist 
Distance: bray 

Dimensions: 2 
Stress:     0.1624969 
Stress type 1, weak ties
Two convergent solutions found after 20 tries
Scaling: centring, PC rotation 
Species: scores missing


ordiplot(Duod.Bray.RI.mds, type = "t")
ordiplot(Duod.Bray.RI.mds, type = "p")


Duod.Classification.RI.map = "Duod.mv.metadata.continuous.RI.txt"
Duod.RI.table=sample_data(read.table(Duod.Classification.RI.map, header=T, sep="\t", row.names=1))
Duod.RI.table

##WARNING: Variables cannot have missing data or n.a, otherwise they'll be counted as categorical variables!!!


Duod.Bray.RI.ef <- envfit(Duod.Bray.RI.mds, Duod.RI.table, permu = 999)
Duod.Bray.RI.ef


#Result
***VECTORS

                           NMDS1    NMDS2     r2 Pr(>r)  
BSS_revised             -0.17629 -0.98434 0.2304  0.039 *
Macrophages              0.20937  0.97784 0.0822  0.353  
Lymphocytes             -0.36364  0.93154 0.1216  0.210  
Neutrophils              0.24641 -0.96917 0.1533  0.152  
Eosinophils              0.08408 -0.99646 0.0823  0.374  
Total_cells             -0.01316 -0.99991 0.0424  0.434  
age_today_dob           -0.00654  0.99998 0.1161  0.233  
RI_impedence_timereflux -0.99999  0.00449 0.0041  0.946  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

plot(Duod.Bray.RI.mds, type = "t", display = "sites")
plot(Duod.Bray.RI.ef, p.max=0.999)


pdf(file = "Duod.Bray.Vectorgram.RI.pdf", width = 12, height = 12)
plot(Duod.Bray.RI.mds, type = "t", display = "sites")
plot(Duod.Bray.RI.ef, p.max=0.999)
dev.off()


#####Calculating univariate distance with adonis for table 2

#Separate by topograph 
#BAL
# use BAL.otu.relative.table
BAL.otu.table = subset_samples(otu.table, Topograph_code==3)

adonis(wUniF.dist ~ age_today_dob*A1, data=data.frame(sample_data(BAL.otu.table)))
sample_data(otu.table)$age_today_dob

sample_data(BAL.otu.table)$age_today_dob

