log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# Load libraries
library(dplyr)
library(biomaRt)

# Define mouse to gene function (for ligands and receptor)
convertMouseGeneList <- function(x){
  # Function obtained from: https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

# Load snakemake inputs.
filtered_degs_table <- snakemake@input[["filtered_degs_table"]]
ligand_receptor_source <- snakemake@params[["ligand_receptor_source"]]
outdir <- snakemake@params[["outdir"]]
sample <- snakemake@params[["sample"]]

# Read filtered degs table.
filtered_degs_table <- read.table(filtered_degs_table, header = TRUE, sep = "\t")

# Read ligands & receptors tables.
ligand_receptor_source_table <- read.table(ligand_receptor_source, header = TRUE, sep = "\t")

human_symbols_table <- convertMouseGeneList(ligand_receptor_source_table$Symbol)
human_symbols <- unlist(lapply(ligand_receptor_source_table$Symbol, function(x) {
  if (x %in% human_symbols_table[,1]){
    gene <- human_symbols_table[human_symbols_table$MGI.symbol == x, "HGNC.symbol"][1]
    print(gene)
  } else {
    gene <- toupper(x)
  }
  return(gene)
}))

ligand_receptor_source_table$Symbol <- human_symbols

# Get ligand receptor gene list.
receptor_list <- filter(ligand_receptor_source_table, Bader_receptor == 1 | FANTOM5_receptor == 1 | CellTalk_receptor == 1)$Symbol
ligand_list <- filter(ligand_receptor_source_table, Bader_ligand == 1 | FANTOM5_ligand == 1 | CellTalk_ligand == 1)$Symbol

# Annotate DEGs table
filtered_degs_table$LR_annotation <- NA
filtered_degs_table$LR_annotation[filtered_degs_table$SYMBOL %in% ligand_list] <- "ligand"
filtered_degs_table$LR_annotation[filtered_degs_table$SYMBOL %in% receptor_list] <- "receptor"
print(head(filtered_degs_table))

# List of ligand receptor belonging to the dataset.
selected_ligands <- subset(filtered_degs_table, LR_annotation == "ligand")$SYMBOL
selected_receptors <- subset(filtered_degs_table, LR_annotation == "receptor")$SYMBOL
merged_vector <- unique(c(selected_ligands, selected_receptors))

# Save files
writeLines(merged_vector, paste0(outdir,"/",sample, "_LR_list.txt"))
write.table(filtered_degs_table, file=paste0(outdir,"/",sample, "_filtered_degs_table_LR_annotated.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
