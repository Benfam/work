setwd 
getwd()
#To load count data and metadata
dat<-read.delim("counts_1.txt", header=T,row.names = 1,sep="\t", comment.char = "!")

info<-read.delim("srp120003_meta_1.txt", header=T,row.names = 1,sep="\t", comment.char = "!")
row.names(dat) -> genes
# Remove prefix from gene

#To link the 2 file to DES
#To construct data into a format for DeSeq2 as test variable


colnames(dat)

#to install DESeq2
BiocManager::install("DESeq2")
library(DESeq2)
head(dat)
head(info)
geneID <- dat$geneid
# Validate sample id for metadata and assay data
all(colnames(dat) == row.names(info))
info$class_id = row.names(info)

info$Type = as.factor(info$Type)

#To construct data into a format for DeSeq2
sa_dds <- DESeqDataSetFromMatrix(countData = dat, 
                                     colData = info, 
                                     design = ~Type)
#normalization of read counts
ddsDE_a = DESeq(sa_dds)
#export normalized read counts
normCounts <-counts(ddsDE_a, normalized = T)
write.csv(normCounts, "normal.counts_1.csv")

#DESeq results
#res_a <-results(ddsDE_a, alpha = 0.01)
res_l <-results(ddsDE_a, contrast = c("Type","DS", "XDR"), alpha = 0.05)

#Export output DESeq results


res_aOrdered2 <- res_l[order(res_l$padj),]
#res_df = res[order(res$log2FoldChange, decreasing = T),]
write.csv(res_aOrdered2, "DESeq3.counts_1.csv")
res_aOrdered2
#summary(res_a)
summary(res_l)

#Visualization
plotMA(ddsDE_a, ylim = c(-3,3))
library(ggplot2)
library(pheatmap)
res_aOrdered2$sig <-ifelse(res_aOrdered2$padj <= 0.01, "yes", "no")

res_deg3 = as.data.frame(res_aOrdered2)
#removing the NA
res_deg3<-na.omit(res_deg3)
### Visualization of DEGs

res_deg3 = as.data.frame(res_aOrdered2)
res_deg3$status = "Not_sig"
res_deg3$status[res_deg3$log2FoldChange>2 & res_deg3$padj < 0.01] = "UP"
res_deg3$status[res_deg3$log2FoldChange < -2 & res_deg3$padj < 0.01] = "DOWN"

## Volcano Plot
library(ggplot2)
ggplot(data = res_deg3, aes(x=log2FoldChange, y=-log10(padj), col=status, label=row.names(res_deg3)))+
  geom_point()+
  theme_minimal()+
  # geom_text_repel()+
  # annotate("text", colour ="black")+
  geom_vline(xintercept = c(-1,1), col="black")+
  geom_hline(yintercept = -log10(0.01), col="black")+
  ggtitle("Resistance vs Susceptible strains ") +
  theme(text = element_text(size = 15))+
  geom_point( size = 3)

#plotMA using log10(baseMean), and log2Foldchange
ggplot(data = res_deg3, aes(x=log10(baseMean), y=log2FoldChange, col=sig, label=row.names(res_deg3)))+
  geom_point()+
  theme_minimal()+
  # geom_text_repel()+
  # annotate("text", colour ="black")+
  geom_vline(xintercept = c(-1,1), col="black")+
  geom_hline(yintercept = -log10(0.01), col="black")+
  ggtitle("Resistance vs Susceptible strains ") +
  theme(text = element_text(size = 15))+
  geom_point( size = 3)


#pheatmap# manipulate the readcounts dataset. Here I filter out DEseq above 0.05 to get only the Significanty expressed genes 
signi<- subset(res_deg3, padj <= 0.01)
# Merge gene for normCounts and DESeq
allsig <- merge(normCounts, signi, by = 0)

write.csv(allsig , "normal_allsig.counts_1.csv")
# to get the normalized read counts
sigCounts =allsig[, 2:13]
#to get the row names and the read counts @ < 0.01
row.names(sigCounts) <- allsig$Row.names

# HeatMap -----------------------------------------------------------------

# data preprocessing 

sig.gene = rownames(res_deg3)[1:15]
expr_heatmap <- normCounts[rownames(normCounts) %in% sig.gene, ]


#Data eady to plot pheatmap showing the expression values
pheatmap(expr_heatmap)

#normalizing the readcounts data using log2 to cater for the outliers
# Here we can see the genes that are highly or lowerly expressed genes within the row as supposed to the raw counts.
#and scale by row to find where the data falls within the distribution/row using standard deviation. Blue is lowerly expressed and red highly expressed while yelow is in between. The heatmap allow us to see how the genes changes between the DS and XRD strains
#pheatmap(log2(expr_heatmap + 1), scale = "row", show_rownames = F, treeheight_row = 0, treeheight_col = 0)

#sigCounts_15 =sigCounts[1:15,]

#pheatmap(sigCounts_15)
#pheatmap(log2(sigCounts_15+ 1), scale = "row", show_rownames = F, treeheight_row = 0, treeheight_col = 0)



#Perform a variance stabilization transformation before a plot
Sample_vst<- vst(sa_dds)

#To perform a variance stabilization transformation before making a plot 
plotPCA(Sample_vst, intgroup = "Type")
#Sample_vst1<- vst(sigCounts)
#plotPCA(sigCounts, intgroup = "Type")

# Significant gene List 

res_deg3_sig = res_deg3[res_deg3$status!="Not_sig",]

# Compute genelist for functional annotation 
library(AnnotationDbi)
library(org.Hs.eg.db)

# Add colunm for Entrez ID to siginificant geen table 
res_deg3_sig$Entrez_IDs=mapIds(org.Hs.eg.db,
                                 keys = row.names(res_deg3_sig),
                                 column = "ENTREZID",
                                 keytype = "SYMBOL",
                                 multiVals = "first")
