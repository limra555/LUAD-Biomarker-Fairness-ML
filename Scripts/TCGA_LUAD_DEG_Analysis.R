# =============================================================================
# TCGA-LUAD: Differential Expression & Pathway Enrichment Analysis
# =============================================================================
# Description : DESeq2-based differential expression analysis of TCGA Lung
#               Adenocarcinoma (LUAD) data, followed by GO and KEGG enrichment.
# Author      : S.Limra Salith
# Date        : 5-5-2025
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load Libraries
# -----------------------------------------------------------------------------

library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)


# -----------------------------------------------------------------------------
# 2. Load Data
# -----------------------------------------------------------------------------

setwd("D:/TCGA-LUAD")

metadata          <- readRDS("metadata.rds")
expression_matrix <- readRDS("expression_matrix.rds")
luad_data         <- readRDS("luad_data.rds")


# -----------------------------------------------------------------------------
# 3. Build DESeqDataSet
# -----------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = expression_matrix,
  colData   = colData(luad_data),
  design    = ~ sample_type
)

# Filter low-count genes (keep genes with total counts >= 10)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set reference level
dds$tissue_type <- factor(dds$tissue_type, levels = c("Normal", "Tumor"))
dds$tissue_type <- relevel(dds$tissue_type, ref = "Normal")


# -----------------------------------------------------------------------------
# 4. Run DESeq2
# -----------------------------------------------------------------------------

dds <- DESeq(dds)


# -----------------------------------------------------------------------------
# 5. Extract & Annotate Results
# -----------------------------------------------------------------------------

res <- results(dds)
res <- as.data.frame(res)
res <- na.omit(res)

# Add Ensembl IDs (strip version numbers)
res$ensembl_gene_id <- sub("\\..*", "", rownames(res))

# Replace zero adjusted p-values to avoid -Inf in log transformation
res$padj[res$padj == 0] <- 1e-300

# Retrieve gene annotations from Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotation <- getBM(
  filters    = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol", "gene_biotype"),
  values     = res$ensembl_gene_id,
  mart       = mart
)

# Merge annotations into results
res <- merge(res, gene_annotation, by = "ensembl_gene_id", all.x = TRUE)


# -----------------------------------------------------------------------------
# 6. Classify Genes by Significance
# -----------------------------------------------------------------------------

res$Significance <- "Not Significant"
res$Significance[res$padj < 0.05 & res$log2FoldChange >  2] <- "Upregulated"
res$Significance[res$padj < 0.05 & res$log2FoldChange < -2] <- "Downregulated"

num_up   <- sum(res$Significance == "Upregulated")
num_down <- sum(res$Significance == "Downregulated")
num_ns   <- sum(res$Significance == "Not Significant")

message("Upregulated genes:     ", num_up)
message("Downregulated genes:   ", num_down)
message("Non-significant genes: ", num_ns)


# -----------------------------------------------------------------------------
# 7. Volcano Plot
# -----------------------------------------------------------------------------

volcano_plot <- ggplot(res,
                       aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  scale_color_manual(values = c(
    "Not Significant" = "grey",
    "Upregulated"     = "red",
    "Downregulated"   = "blue"
  )) +
  
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  annotate("text",
           x     = max(res$log2FoldChange, na.rm = TRUE),
           y     = max(-log10(res$padj),   na.rm = TRUE),
           label = paste("Upregulated:", num_up),
           color = "red", hjust = 1
  ) +
  
  annotate("text",
           x     = min(res$log2FoldChange, na.rm = TRUE),
           y     = max(-log10(res$padj),   na.rm = TRUE),
           label = paste("Downregulated:", num_down),
           color = "blue", hjust = 0
  ) +
  
  labs(
    title = "Volcano Plot of Differentially Expressed Genes (LUAD)",
    x     = "Log2 Fold Change",
    y     = "-Log10 Adjusted P-value"
  ) +
  
  theme_minimal()

print(volcano_plot)

ggsave("VolcanoPlot_LUAD.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)


# -----------------------------------------------------------------------------
# 8. Filter DEGs
# -----------------------------------------------------------------------------

# All significant DEGs (any biotype)
deg_filtered <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 2)
message("Total significant DEGs: ", nrow(deg_filtered))

# Protein-coding DEGs only
deg_protein_coding <- subset(deg_filtered, gene_biotype == "protein_coding")
message("Protein-coding DEGs: ", nrow(deg_protein_coding))

# Split by direction
up_protein   <- subset(deg_protein_coding, log2FoldChange >= 2)
down_protein <- subset(deg_protein_coding, log2FoldChange <= -2)

message("Protein-coding Upregulated:   ", nrow(up_protein))
message("Protein-coding Downregulated: ", nrow(down_protein))


# -----------------------------------------------------------------------------
# 9. Save DEG Tables
# -----------------------------------------------------------------------------

write.csv(res,               "All_DEGs_LUAD.csv",                    row.names = FALSE)
write.csv(deg_filtered,      "Significant_DEGs_LUAD.csv",            row.names = FALSE)
write.csv(deg_protein_coding,"ProteinCoding_DEGs_LUAD.csv",          row.names = FALSE)
write.csv(up_protein,        "Upregulated_ProteinCoding_DEGs.csv",   row.names = FALSE)
write.csv(down_protein,      "Downregulated_ProteinCoding_DEGs.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 10. Pathway Enrichment — Prepare Gene Lists
# -----------------------------------------------------------------------------

options(timeout = 300)   # Increase timeout for KEGG downloads

up_genes   <- unique(up_protein$hgnc_symbol[!is.na(up_protein$hgnc_symbol)])
down_genes <- unique(down_protein$hgnc_symbol[!is.na(down_protein$hgnc_symbol)])

# Convert gene symbols to Entrez IDs
up_entrez <- bitr(up_genes,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
down_entrez <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


# -----------------------------------------------------------------------------
# 11. GO Enrichment (Biological Process)
# -----------------------------------------------------------------------------

go_up <- enrichGO(
  gene          = up_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

go_down <- enrichGO(
  gene          = down_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)


# -----------------------------------------------------------------------------
# 12. KEGG Enrichment
# -----------------------------------------------------------------------------

kegg_up <- enrichKEGG(
  gene         = up_entrez$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)

kegg_down <- enrichKEGG(
  gene         = down_entrez$ENTREZID,
  organism     = "hsa",
  pvalueCutoff = 0.05
)


# -----------------------------------------------------------------------------
# 13. Save Enrichment Tables
# -----------------------------------------------------------------------------

write.csv(as.data.frame(go_up),   "GO_BP_Upregulated_LUAD.csv",   row.names = FALSE)
write.csv(as.data.frame(go_down), "GO_BP_Downregulated_LUAD.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_up), "KEGG_Upregulated_LUAD.csv",    row.names = FALSE)
write.csv(as.data.frame(kegg_down),"KEGG_Downregulated_LUAD.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 14. Plot Enrichment Results
# -----------------------------------------------------------------------------

# GO barplots
barplot(go_up,
        showCategory = 10,
        title = "GO Biological Process – Upregulated Genes"
)
ggsave("GO_BP_Upregulated.png", width = 8, height = 6, dpi = 300)

barplot(go_down,
        showCategory = 10,
        title = "GO Biological Process – Downregulated Genes"
)
ggsave("GO_BP_Downregulated.png", width = 8, height = 6, dpi = 300)

# KEGG dotplots
dotplot(kegg_up,
        showCategory = 10,
        title = "KEGG Pathways – Upregulated Genes"
)
ggsave("KEGG_Upregulated.png", width = 8, height = 6, dpi = 300)

dotplot(kegg_down,
        showCategory = 10,
        title = "KEGG Pathways – Downregulated Genes"
)
ggsave("KEGG_Downregulated.png", width = 8, height = 6, dpi = 300)

# =============================================================================
# End of Script
# =============================================================================