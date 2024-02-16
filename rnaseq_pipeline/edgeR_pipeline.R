
library(limma)
library(edgeR)
#library(EDASeq)
#library(RUVSeq)
#library(ffpe)
library(RColorBrewer)
library(tidyverse)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(pheatmap)
library(EnhancedVolcano)

#sshfs yumei@nrnb-login.ucsd.edu:/cellar/users/yumei/Documents/nlp_mRNA/experimental_RNAseq_analysis/mouse_lji_110623_analysis/rsem_output /Users/lb/Documents/remote_access


load("/Users/lb/Documents/remote_access/RSEM_Quant.genes.counts.RData")

#loading the data:========================

all_samples <- counts
colnames(all_samples) <- gsub("_rsem.*", "", colnames(all_samples))


#specifying the conditions================
group <- c("a", "a", "a", "a", "a", "a", "a", "b", "b", "b", 'b', "b", "b", "b", "c", "c", "c", "c", "c", "c")

fraction <- factor(c("Cytosol","Cytosol", "Cytosol", "Cytosol", "Cytosol", "Cytosol", "Cytosol", "Synapse","Synapse","Synapse","Synapse", "Synapse","Synapse","Synapse", 'Lysate', 'Lysate', 'Lysate', 'Lysate', 'Lysate', 'Lysate'))

condition_colors <- c('gray', "gray", "gray", "gray", "gray", "gray", "gray", "red", "red", "red", "red", "red", "red", "red", "blue", "blue",  "blue",  "blue", "blue",  "blue")

data.frame(Sample=colnames(all_samples), fraction)

group <-paste0(fraction)

design <- model.matrix(~0+group)

cpm_log <- cpm(all_samples, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)

#filtering:
expr_cutoff <- 1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data_clean <- all_samples[median_log2_cpm > expr_cutoff, ]
#13456

cpm_log <- cpm(data_clean, log = TRUE)

corr_df <- cor(cpm_log)

annotation_row = data.frame(Fraction=fraction)
rownames(annotation_row) = rownames(corr_df)
pheatmap(corr_df, annotation_row=annotation_row)

#plot the top 500 most variable genes heatmap across conditions
var_genes <- apply(cpm_log, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
highly_variable_lcpm <- cpm_log[select_var,]
var_corr_df <- cor(highly_variable_lcpm)
annotation_row = data.frame(Fraction=fraction)
rownames(annotation_row) = rownames(var_corr_df)
pheatmap(var_corr_df, annotation_row=annotation_row)



# # 
pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = 16, col=condition_colors, cex= 3, xlab = "PC1", ylab = "PC2")
#text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
summary(pca)

unique_condition_colors <- unique(condition_colors)

legend(x='bottomleft', legend=c("Cytosol", "Synapse", "Total Lysate"), 
       fill = unique_condition_colors)
)

# Perform PCA and store the summary
pca <- prcomp(t(cpm_log), scale. = TRUE)
pca_summary <- summary(pca)

# Calculate the percentage of variance explained by PC1 and PC2
pc1_variance <- pca_summary$sdev[1]^2 / sum(pca_summary$sdev^2) * 100
pc2_variance <- pca_summary$sdev[2]^2 / sum(pca_summary$sdev^2) * 100

# Create the plot, including the percentage of variance explained by PC1 and PC2 in the axis labels
plot(pca$x[, 1], pca$x[, 2], pch = 16, col=condition_colors, cex= 3, 
     xlab = paste0("PC1: ", round(pc1_variance, 2), "% variance explained"),
     ylab = paste0("PC2: ", round(pc2_variance, 2), "% variance explained"))


y <- DGEList(counts=data_clean, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

my.contrasts <- makeContrasts(SynapsevsCytosol=groupSynapse-groupCytosol, SynapsevsLysate=groupSynapse-groupLysate, levels=design)

#compare Synapse vs Cytosol:================
qlf.SynapsevsCytosol <- glmQLFTest(fit, contrast=my.contrasts[, "SynapsevsCytosol"])
topTags(qlf.SynapsevsCytosol)
results_edgeR <- topTags(qlf.SynapsevsCytosol, n = nrow(data_clean), sort.by = "none")
sum(results_edgeR$table$FDR < .05)
plotSmear(qlf.SynapsevsCytosol, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .05])
#11454 DEGs at FDR <0.05
library("org.Mm.eg.db")
results.df <- as.data.frame(results_edgeR)
head(results.df)
results.df$symbol <- mapIds(org.Mm.eg.db, keys=rownames(results.df), keytype= "ENSEMBL", column='SYMBOL')
write.csv(results.df, '/Users/lb/Documents/mouse_052323_newseq_analysis/output_files/052323_newseq_SynapsevsCytosol.csv')
# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  results.df$logFC < -1 & results.df$FDR<0.05, 'red',
  ifelse(results.df$logFC > 1 & results.df$FDR<0.05, 'darkgreen',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'darkgreen'] <- 'high'
names(keyvals)[keyvals == 'gray'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'

# Now, pass this column to the lab argument of the EnhancedVolcano function
EnhancedVolcano(results.df, lab =results.df$symbol,
                selectLab = c('Nlgn2', 'Scn8a', 'Cacna1i', 'Nrxn2', 'Grm3', 'Cacna1b', 'Nrxn3', 'Celsr2','Grm1', 'Grm2', 'Grm4', 'Grm5'),
                x = 'logFC', 
                y = 'FDR', 
                pointSize = 2, 
                labSize = 5, 
                FCcutoff = 1, 
                pCutoff = 0.05, 
                colCustom = keyvals, 
                colAlpha = 1)


#GO Enrichment:
deg <- results.df[results.df$FDR < .05 & results.df$logFC>0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

deg <- results.df[results.df$FDR < .05 & results.df$logFC<0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)



#Compare Synapse vs. Lysate:============================================================
qlf.SynapsevsLysate <- glmQLFTest(fit, contrast=my.contrasts[, "SynapsevsLysate"])
topTags(qlf.SynapsevsLysate)
results_edgeR <- topTags(qlf.SynapsevsLysate, n = nrow(data_clean), sort.by = "none")
sum(results_edgeR$table$FDR < .05)
plotSmear(qlf.SynapsevsLysate, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .05])
#10088 DEGs at FDR <0.05; 
library("org.Mm.eg.db")
results.df <- as.data.frame(results_edgeR)
head(results.df)
results.df$symbol <- mapIds(org.Mm.eg.db, keys=rownames(results.df), keytype= "ENSEMBL", column='SYMBOL')

write.csv(results.df, '/Users/lb/Documents/mouse_052323_newseq_analysis/output_files/052323_newseq_SynapsevsLysate.csv')

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  results.df$logFC < -1 & results.df$FDR<0.05, 'red',
  ifelse(results.df$logFC > 1 & results.df$FDR<0.05, 'darkgreen',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'darkgreen'] <- 'high'
names(keyvals)[keyvals == 'gray'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'


#EnhancedVolcano(results.df, lab =results.df$symbol, x='logFC', y='PValue', pointSize=2, labSize=6, FCcutoff=1, pCutoff=0.0274, ylim = c(0, -log10(10e-12)))
EnhancedVolcano(results.df, lab =results.df$symbol, x='logFC', y='FDR', pointSize=2, labSize=6, FCcutoff=1, pCutoff=0.05, colCustom=keyvals, colAlpha=1, selectLab = c('Nlgn2', 'Scn8a', 'Cacna1i', 'Nrxn2', 'Grm3', 'Cacna1b', 'Nrxn3', 'Celsr2', 'Grm1', 'Grm2', 'Grm4', 'Grm5'))
#pCutoff=0.001 is equal to FDR = 0.15


#GO Enrichment:
deg <- results.df[results.df$FDR < .05 & results.df$logFC>0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

deg <- results.df[results.df$FDR < .05 & results.df$logFC<0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

#==isoforms===================================================================================================================================
load("/Users/lb/Documents/remote_access/RSEM_Quant.isoforms.counts.RData")
all_samples <- counts
colnames(all_samples) <- gsub("_rsem.*", "", colnames(all_samples))


#specifying the conditions================
group <- c("a", "a", "a", "a", "a", "a", "a", "b", "b", "b", 'b', "b", "b", "b", "c", "c", "c", "c", "c", "c")

fraction <- factor(c("Cytosol","Cytosol", "Cytosol", "Cytosol", "Cytosol", "Cytosol", "Cytosol", "Synapse","Synapse","Synapse","Synapse", "Synapse","Synapse","Synapse", 'Lysate', 'Lysate', 'Lysate', 'Lysate', 'Lysate', 'Lysate'))

condition_colors <- c('gray', "gray", "gray", "gray", "gray", "gray", "gray", "red", "red", "red", "red", "red", "red", "red", "blue", "blue",  "blue",  "blue", "blue",  "blue")

data.frame(Sample=colnames(all_samples), fraction)

group <-paste0(fraction)

design <- model.matrix(~0+group)

cpm_log <- cpm(all_samples, log = TRUE)
median_log2_cpm <- apply(cpm_log, 1, median)
hist(median_log2_cpm)

#filtering:
expr_cutoff <- 1
abline(v = expr_cutoff, col = "red", lwd = 3)
sum(median_log2_cpm > expr_cutoff)
data_clean <- all_samples[median_log2_cpm > expr_cutoff, ]
#13456

cpm_log <- cpm(data_clean, log = TRUE)

corr_df <- cor(cpm_log)

annotation_row = data.frame(Fraction=fraction)
rownames(annotation_row) = rownames(corr_df)
pheatmap(corr_df, annotation_row=annotation_row)

#plot the top 500 most variable genes heatmap across conditions
var_genes <- apply(cpm_log, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
highly_variable_lcpm <- cpm_log[select_var,]
var_corr_df <- cor(highly_variable_lcpm)
annotation_row = data.frame(Fraction=fraction)
rownames(annotation_row) = rownames(var_corr_df)
pheatmap(var_corr_df, annotation_row=annotation_row)



# # 
pca <- prcomp(t(cpm_log), scale. = TRUE)
plot(pca$x[, 1], pca$x[, 2], pch = 16, col=condition_colors, cex= 3, xlab = "PC1", ylab = "PC2")
#text(pca$x[, 1], pca$x[, 2], labels = colnames(cpm_log))
summary(pca)

unique_condition_colors <- unique(condition_colors)

legend(x='bottomleft', legend=c("Cytosol", "Synapse", "Total Lysate"),fill = unique_condition_colors)

# Perform PCA and store the summary
pca <- prcomp(t(cpm_log), scale. = TRUE)
pca_summary <- summary(pca)

# Calculate the percentage of variance explained by PC1 and PC2
pc1_variance <- pca_summary$sdev[1]^2 / sum(pca_summary$sdev^2) * 100
pc2_variance <- pca_summary$sdev[2]^2 / sum(pca_summary$sdev^2) * 100

# Create the plot, including the percentage of variance explained by PC1 and PC2 in the axis labels
plot(pca$x[, 1], pca$x[, 2], pch = 16, col=condition_colors, cex= 3, 
     xlab = paste0("PC1: ", round(pc1_variance, 2), "% variance explained"),
     ylab = paste0("PC2: ", round(pc2_variance, 2), "% variance explained"))


y <- DGEList(counts=data_clean, group=group)
y <- calcNormFactors(y)
y <- estimateDisp(y,design)
fit <- glmQLFit(y, design)

my.contrasts <- makeContrasts(SynapsevsCytosol=groupSynapse-groupCytosol, SynapsevsLysate=groupSynapse-groupLysate, levels=design)

#compare Synapse vs Cytosol:================
qlf.SynapsevsCytosol <- glmQLFTest(fit, contrast=my.contrasts[, "SynapsevsCytosol"])
topTags(qlf.SynapsevsCytosol)
results_edgeR <- topTags(qlf.SynapsevsCytosol, n = nrow(data_clean), sort.by = "none")
sum(results_edgeR$table$FDR < .05)
plotSmear(qlf.SynapsevsCytosol, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .05])
#11454 DEGs at FDR <0.05
library("org.Mm.eg.db")
# 10041 DEGs at FDR <0.05
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

results.df <- as.data.frame(results_edgeR)

genes <- getBM(attributes = c('ensembl_transcript_id', 'external_gene_name'), 
               filters = 'ensembl_transcript_id', 
               values = rownames(results.df), 
               mart = mart)

results.df <- merge(results.df, genes, by.x = "row.names", by.y = "ensembl_transcript_id", all.x = TRUE)
row.names(results.df) <- results.df$Row.names
results.df$Row.names <- NULL

head(results.df)
write.csv(results.df, '/Users/lb/Documents/mouse_052323_newseq_analysis/output_files/052323_newseq_isoforms_SynapsevsCytosol.csv')


# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  results.df$logFC < -1 & results.df$FDR<0.05, 'red',
  ifelse(results.df$logFC > 1 & results.df$FDR<0.05, 'darkgreen',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'darkgreen'] <- 'high'
names(keyvals)[keyvals == 'gray'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'

# Now, pass this column to the lab argument of the EnhancedVolcano function
EnhancedVolcano(results.df, lab =results.df$symbol,
                selectLab = c('Nlgn2', 'Scn8a', 'Cacna1i', 'Nrxn2', 'Grm3', 'Cacna1b', 'Nrxn3', 'Celsr2','Grm1', 'Grm2', 'Grm4', 'Grm5'),
                x = 'logFC', 
                y = 'FDR', 
                pointSize = 2, 
                labSize = 5, 
                FCcutoff = 1, 
                pCutoff = 0.05, 
                colCustom = keyvals, 
                colAlpha = 1)


#GO Enrichment:
deg <- results.df[results.df$FDR < .05 & results.df$logFC>0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

deg <- results.df[results.df$FDR < .05 & results.df$logFC<0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)



#Compare Synapse vs. Lysate:============================================================
qlf.SynapsevsLysate <- glmQLFTest(fit, contrast=my.contrasts[, "SynapsevsLysate"])
topTags(qlf.SynapsevsLysate)
results_edgeR <- topTags(qlf.SynapsevsLysate, n = nrow(data_clean), sort.by = "none")
sum(results_edgeR$table$FDR < .05)
plotSmear(qlf.SynapsevsLysate, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < .05])
#10088 DEGs at FDR <0.05; 
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

results.df <- as.data.frame(results_edgeR)

genes <- getBM(attributes = c('ensembl_transcript_id', 'external_gene_name'), 
               filters = 'ensembl_transcript_id', 
               values = rownames(results.df), 
               mart = mart)

results.df <- merge(results.df, genes, by.x = "row.names", by.y = "ensembl_transcript_id", all.x = TRUE)
row.names(results.df) <- results.df$Row.names
results.df$Row.names <- NULL

head(results.df)
write.csv(results.df, '/Users/lb/Documents/mouse_052323_newseq_analysis/output_files/052323_newseq_isoforms_SynapsevsLysate.csv')

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  results.df$logFC < -1 & results.df$FDR<0.05, 'red',
  ifelse(results.df$logFC > 1 & results.df$FDR<0.05, 'darkgreen',
         'gray'))
keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'darkgreen'] <- 'high'
names(keyvals)[keyvals == 'gray'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'


#EnhancedVolcano(results.df, lab =results.df$symbol, x='logFC', y='PValue', pointSize=2, labSize=6, FCcutoff=1, pCutoff=0.0274, ylim = c(0, -log10(10e-12)))
EnhancedVolcano(results.df, lab =results.df$symbol, x='logFC', y='FDR', pointSize=2, labSize=6, FCcutoff=1, pCutoff=0.05, colCustom=keyvals, colAlpha=1, selectLab = c('Nlgn2', 'Scn8a', 'Cacna1i', 'Nrxn2', 'Grm3', 'Cacna1b', 'Nrxn3', 'Celsr2', 'Grm1', 'Grm2', 'Grm4', 'Grm5'))
#pCutoff=0.001 is equal to FDR = 0.15


#GO Enrichment:
deg <- results.df[results.df$FDR < .05 & results.df$logFC>0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)

deg <- results.df[results.df$FDR < .05 & results.df$logFC<0,]
ego2 <- enrichGO(gene         = deg$symbol,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(summary(ego2))

dotplot(ego2, showCategory=15)



