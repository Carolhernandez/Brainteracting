log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Load libraries
library(dplyr)

# Load snakemake inputs.
graph <- snakemake@input[["graph"]]
filtered_degs_tables_ann <- snakemake@input[["filtered_degs_tables_ann"]]
PPI_score_threshold <- snakemake@params[["PPI_score_threshold"]]
outdir <- snakemake@params[["outdir"]]
samples <- snakemake@params[["samples"]]

# Read DEGs table and set the names.
degs_table_list <- lapply(1:2, function(x) {
  degs_table <- read.table(filtered_degs_tables_ann[x], header = TRUE, sep = "\t")
  degs_table_filt <- degs_table[!(is.na(degs_table$LR_annotation)),c("SYMBOL", "LR_annotation")]
  degs_table_filt$tissue <- samples[x]
  return(degs_table_filt)
})
names(degs_table_list) <- samples

# Read the PPI network.
graph <- read.table(graph, sep = "\t", header = TRUE)


# Create node annotation file.
merge <- merge(degs_table_list[[1]], degs_table_list[[2]], by = "SYMBOL", all = TRUE, suffixes = paste0(".", samples)) # DEGs dataframe are merged
row.names(merge) <- merge$SYMBOL
LR_ann <- sapply(merge$SYMBOL, function(x) {ifelse(is.na(merge[x,2]), merge[x,4], merge[x,2])}) # Ligand-receptor vector is obtained
tissue_type <- rep("both", dim(merge)[1]) # Tissue vector is obtained
tissue_type[is.na(merge[,2])] <- samples[2] # Tissue vector is obtained
tissue_type[is.na(merge[,4])] <- samples[1] # Tissue vector is obtained
merge$LR_annotation <- LR_ann
merge$tissue_type <- tissue_type
node_annotation <- merge[,c(1,6,7)] # Select new tables.


# Tissue check.
# Get inputs
genes_sample1 <- degs_table_list[[1]]$SYMBOL
genes_sample2 <- degs_table_list[[2]]$SYMBOL
both_genes <- subset(node_annotation, tissue_type == "both")$SYMBOL
bool_vec_tissue <- rep(NA, dim(graph)[1])

#Apply filter.
for (i in 1:length(bool_vec_tissue)){
  proteins <- graph[i, c("protein1", "protein2")]
  value <- FALSE
  if (proteins[1] %in% genes_sample1 & proteins[2] %in% genes_sample2) {value <- TRUE}
  if (proteins[2] %in% genes_sample1 & proteins[1] %in% genes_sample2) {value <- TRUE}
  if (proteins[1] %in% both_genes | proteins[2] %in% both_genes) {value <- TRUE}
  bool_vec_tissue[i] <- value
}
tissue_filtered_graph <- graph[bool_vec_tissue,]


# LR check.
# Get inputs
receptors <- subset(node_annotation, LR_annotation == "receptor")$SYMBOL
ligands <- subset(node_annotation, LR_annotation == "ligand")$SYMBOL
bool_vec_LR <- rep(NA, dim(tissue_filtered_graph)[1])

#Apply filter.
for (i in 1:length(bool_vec_LR)){
  proteins <- tissue_filtered_graph[i, c("protein1", "protein2")]
  value <- FALSE
  if (proteins[1] %in% receptors & proteins[2] %in% ligands) {value <- TRUE}
  if (proteins[1] %in% ligands & proteins[2] %in% receptors) {value <- TRUE}
  bool_vec_LR[i] <- value
}
complete_filtered_graph <- tissue_filtered_graph[bool_vec_LR,]

# Filter PPI by minimum threshold.
PPI_score_filtered_graph <- filter(complete_filtered_graph, combined_score >= PPI_score_threshold)

# Filter previous annotation file by filtered PPI nodes.
filtered_node_annotation <- node_annotation[node_annotation$SYMBOL %in% unique(c(PPI_score_filtered_graph$protein1, PPI_score_filtered_graph$protein1)),]


# Write outputs.
write.table(PPI_score_filtered_graph, file = paste0(outdir, "/PPI_network_filtered.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filtered_node_annotation, file = paste0(outdir, "/node_annotation_file.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
