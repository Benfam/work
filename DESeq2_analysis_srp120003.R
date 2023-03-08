
setwd 
getwd()
#To load count data and metadata
SRA_counts<-read.delim("counts_1.txt", header=T,row.names = 1,sep="\t", comment.char = "!")

pheno<-read.delim("srp120003_meta_1.txt", header=T,row.names = 1,sep="\t", comment.char = "!")
row.names(SRA_counts) -> genes
# Remove prefix from gene



colnames(SRA_counts)

#to install DESeq2
BiocManager::install("DESeq2")
library(DESeq2)
head(SRA_counts)
head(pheno)
geneID <- SRA_counts$geneid
# Validate sample id for metadata and assay data
all(colnames(SRA_counts) == row.names(pheno))
pheno$class_id = row.names(pheno)

pheno$Type = as.factor(pheno$Type)


#SRA_counts <-SRA_counts[, as.character(pheno$Sample.Name)] 
#colnames(SRA_counts) == pheno$Sample.Name, na.rm = FALSE)

#To construct data into a format for DeSeq2
sample_dds <- DESeqDataSetFromMatrix(countData = SRA_counts, 
                                     colData = pheno, 
                                     design = ~Type)

dds = DESeq(sample_dds)

res <-results(dds, contrast = c("Type","DS", "XDR"))
#res_t <-results(dds, alpha = 0.01)
#summary(res)

res_df = res[order(res$log2FoldChange, decreasing = T),]

# Visulise DEGs
res_deg = as.data.frame(res_df)
res_deg$status = "Not_sig"
res_deg$status[res_deg$log2FoldChange>2 & res_deg$padj < 0.01] = "UP"
res_deg$status[res_deg$log2FoldChange < -2 & res_deg$padj < 0.01] = "DOWN"

library(ggplot2)
ggplot(data = res_deg, aes(x=log2FoldChange, y=-log10(padj), col=status, label=row.names(res_deg)))+
  geom_point()+
  theme_minimal()+
  # geom_text_repel()+
  # annotate("text", colour ="black")+
  geom_vline(xintercept = c(-1,1), col="black")+
  geom_hline(yintercept = -log10(0.01), col="black")+
  ggtitle("Resistance vs Susceptible strains ") +
  theme(text = element_text(size = 20))+
  geom_point( size = 3)


#Perform a variance stabilization transformation before a plot
Sample_vst<- vst(sample_dds)
#To perform a variance stabilization transformation before making a plot 
plotPCA(Sample_vst, intgroup = "Type")
