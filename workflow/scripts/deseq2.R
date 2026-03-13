# DESeq2 differential expression analysis
# Called by Snakemake rule; inputs/outputs/params passed via snakemake object

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")

library(DESeq2)
library(openxlsx)

htseq_dir       <- snakemake@params$htseq_dir
samples_tsv     <- snakemake@input$samples
condition_col   <- snakemake@params$condition_col
reference_level <- snakemake@params$reference_level
out_results     <- snakemake@output$results
out_rdata       <- snakemake@output$rdata

# Load sample sheet
coldata <- read.table(samples_tsv, header = TRUE, sep = "\t", row.names = "sample")

# Read HTSeq count files
count_files <- file.path(htseq_dir, paste0(rownames(coldata), ".tsv"))
names(count_files) <- rownames(coldata)

counts_list <- lapply(count_files, function(f) {
  d <- read.table(f, header = FALSE, row.names = 1)
  d[!grepl("^__", rownames(d)), , drop = FALSE]
})
count_matrix <- do.call(cbind, counts_list)
colnames(count_matrix) <- rownames(coldata)

# Set reference level
coldata[[condition_col]] <- relevel(
  factor(coldata[[condition_col]]),
  ref = reference_level
)

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = coldata,
  design    = as.formula(paste("~", condition_col))
)
dds <- DESeq(dds)
save(dds, file = out_rdata)

# Extract and export results
res <- results(dds, alpha = 0.05)
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]

dir.create(dirname(out_results), showWarnings = FALSE, recursive = TRUE)
write.xlsx(res_df, out_results, rowNames = FALSE)
message("DESeq2 complete. Results written to: ", out_results)
