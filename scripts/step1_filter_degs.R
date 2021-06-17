log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Load libraries
library(dplyr)
library(biomaRt)

# Mouse to human function is defined.
convertMouseGeneList <- function(x){
  # Function obtained from: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

# Load snakemake inputs.
degs_table <- snakemake@input[["degs_table"]]
filter_table <- snakemake@params[["filter_table"]]
outdir <- snakemake@params[["outdir"]]
sample <- snakemake@params[["sample"]]
species <- snakemake@params[["species"]]

# DEGs table is loaded.
degs_table <- read.table(degs_table, header = TRUE, sep = ",")

# Samples species check.
if (species == "mouse"){
  human_symbols_table <- convertMouseGeneList(degs_table$SYMBOL) # Get table.

  # Apply for each genes where (1) mouse to human (2) mouse to upper font.
  human_symbols <- unlist(lapply(degs_table$SYMBOL, function(x) {
    if (x %in% human_symbols_table[,1]){
      gene <- human_symbols_table[human_symbols_table$MGI.symbol == x, "HGNC.symbol"][1]
    } else {
      gene <- toupper(x)
    }
    return(gene)
  }))

  # Update symbols.
  degs_table$SYMBOL <- human_symbols
}

# Get filtering values.
filter_table <- read.table(filter_table, header = TRUE, sep = "\t")
FDR_threshold <- filter_table[filter_table$sample == sample, "FDR"]
logFC_threshold <- filter_table[filter_table$sample == sample, "logFC"]

# Filter DEGs tables.
filtered_degs_table <- filter(degs_table, padj <= FDR_threshold & log2FoldChange >= logFC_threshold & log2FoldChange >= 0)

# Save outputs.
write.table(filtered_degs_table, file=paste0(outdir,"/",sample, "_filtered_degs_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
