# Libraries

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")
BiocManager::install("pathview")

# Load the count data
count_data <- read.csv("C:/Users/ignac/Desktop/R/Curso/Basics/count_matrix.csv", sep = ";", header = TRUE, row.names = 1)
View(count_data)
colnames(count_data)
head(count_data)

#Load the sample data 
sample_info <- read.csv("C:/Users/ignac/Desktop/R/Curso/Basics/design2.csv", sep = ";")
View(sample_info)
colnames(sample_info)
head(sample_info)

#Set factors
sample_info$Treatment <- factor(sample_info$Treatment)
sample_info$Sequencing <- factor(sample_info$Sequencing)

# Create a deseq object and import the count data and the sample information
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = sample_info, design = ~ Sequencing + Treatment)
View(dds)

#Set the reference for the treatment factor
dds$Treatment <- factor(dds$Treatment, levels = c("untreated", "treated"))

#Filter the genes
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

#Perform statistical test for identifying DEGs
dds <- DESeq(dds)
deseq_result <- results(dds)
deseq_result

#Convert in R dataframe
deseq_result <- as.data.frame(deseq_result)
class(deseq_result)
View(deseq_result)

#Order by increasing p-value
deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
View(deseq_result_ordered)

#See one specific gene
deseq_result_ordered["FBgn0261552",]

#Filter by genes that have less the 0.05 in padj and < -1 and > 1
filtered <- deseq_result %>% filter(deseq_result$padj < 0.05 & abs(deseq_result$log2FoldChange) > 1)
View(filtered)

#Save results
write.csv(deseq_result, "de_result_all.csv")
write.csv(filtered, "de_result_filtered.csv")
normalized_counts <- counts(dds, normalized = TRUE)
View(normalized_counts)
write.csv(normalized_counts, "normalized_counts.csv")

#Visualization
#Dispersion plot
plotDispEsts(dds)

#PCA
#variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup=c("Sequencing", "Treatment"))

#Heatmap
#Distance matrix heatmap
#Generate the distance matrix
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
View(sampleDistMatrix)
#Generate heatmap
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist, col = colors)

#Heatmap of log transformed normalized counts
top_hits <- deseq_result[order(deseq_result$padj),][1:10,]
top_hits <- row.names(top_hits)
top_hits

rld <- rlog(dds, blind = FALSE)
pheatmap(assay(rld)[top_hits,])

#Heatmap of Z score -> Top 10 genes
cal_z_score <- function(x) {(x - mean(x))/ sd(x)}
zscore_all <- t(apply(normalized_counts,1,cal_z_score))
zscore_subset <- zscore_all[top_hits,]
pheatmap(zscore_subset)


#MAplot
plotMA(dds, ylim=c(-2,2))

#Remove the noise
resLFC <- lfcShrink(dds, coef = "Treatment_treated_vs_untreated", type = "apeglm")
plotMA(resLFC, ylim=c(-2,2))

#Volcano Plot
#resLFC to data frame
resLFC <- as.data.frame(resLFC)
#Label the genes

resLFC$diffexpressed <- "NO"
resLFC$diffexpressed[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05] <- "DOWN"
resLFC$delabel <- NA

ggplot(data = resLFC, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed, label = delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  theme(text = element_text(size =20))

#KEGG pathway

library(tibble)
df <- tibble::rownames_to_column(filtered, "VALUE")
colnames(df)[1] = "GENEID"

filtered$ENTREZ = mapIds(org.Mm.eg.db,
                         key = row.names(filtered),
                         column = "ENTREZID",
                         multiVals = "first")

foldchange = filtered$log2FoldChange
names(foldchange) = row.names(filtered)
head(foldchange)
View(foldchange)

data("go.sets.mm")
data("go.subs.mm")

gobpsets = go.sets.mm[go.subs.mm$BP]
gobpres = gage(exprs = foldchange, gsets = gobpsets, same.dir = TRUE)
View(gobpres)
