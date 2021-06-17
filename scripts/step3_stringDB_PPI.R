log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Load libraries
library(STRINGdb)
library(biomaRt)
library(dplyr)

# Load snakemake inputs.
input_genes_lists <- snakemake@input[["LR_lists"]]
mode <- snakemake@params[["mode"]]
stringdb_path <- snakemake@params[["stringDB"]]
outdir <- snakemake@params[["outdir"]]


# 1. RL lists to unique vector.
input_genes <- unique(unlist(lapply(input_genes_lists, function(x) {scan(what = "character", file = x)})))

# 2. Get PPI network.
# User recommendation not to use mode online. Instead, download he PPI interactions file from string database and continue with offline mode.
if (mode == "online"){
  example_df <- data.frame("logFC" = rep(1, 341), "FDR" = rep(1, 341), gene = input_genes)
  string_db <- STRINGdb$new( version="11", species=9606, score_threshold=200, input_directory="")
  example1_mapped <- string_db$map(example_df, "gene", removeUnmappedRows = TRUE)
  interactions <- string_db$get_interactions(example1_mapped[,4])
  protein1 <- sapply(interactions$from, function(x) {example1_mapped[example1_mapped$STRING_id == x, "gene"]})
  protein2 <- sapply(interactions$to, function(x) {example1_mapped[example1_mapped$STRING_id == x, "gene"]})
  interactions$from <- protein1
  interactions$to <- protein2
  interactions_filtered <- interactions %>% distinct()
}

if (mode == "offline"){
  # Read offline stringDB
  string_db_offline <- read.table(file = stringdb_path, sep = " ", header = TRUE)

  #Get biomart
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

  # Get biomart translator.
  translator <- getBM(attributes = c("ensembl_peptide_id", "hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values = input_genes,
                      mart = ensembl)

  # Extract all protein ID for input genes and paste stringDB identifier.
  input_genes_prot_ID <- paste(9606, unique(translator$ensembl_peptide_id), sep = ".")

  # Filter stringDB by input genes.
  filtered_string_db <- string_db_offline[(string_db_offline$protein1 %in% input_genes_prot_ID & string_db_offline$protein2 %in% input_genes_prot_ID),]

  # Recover HUGO IDs
  recovered_IDs_1 <- sapply(strsplit(filtered_string_db$protein1, "\\."), `[[`, 2)
  recovered_IDs_2 <- sapply(strsplit(filtered_string_db$protein2, "\\."), `[[`, 2)
  hugo_protein1 <- unlist(lapply(recovered_IDs_1, function(x) {return(subset(translator, ensembl_peptide_id == x)$hgnc_symbol)}))
  hugo_protein2 <- unlist(lapply(recovered_IDs_2, function(x) {return(subset(translator, ensembl_peptide_id == x)$hgnc_symbol)}))
  filtered_string_db$protein1 <- hugo_protein1
  filtered_string_db$protein2 <- hugo_protein2
  filtered_string_db_unique <- filtered_string_db %>% distinct()

}

# 3. Write PPI table.
write.table(filtered_string_db_unique, file = paste0(outdir, "/PPI_network.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
