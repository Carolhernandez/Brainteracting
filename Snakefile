import pandas as pd
from snakemake.utils import min_version
##### set minimum snakemake version #####
min_version("5.10.0")

configfile: "config.yaml"

OUTDIR = config["outdir"]
LOGDIR = config["logdir"]


samples = pd.read_csv(config["samples"], sep="\t", comment="#")
samples = samples.set_index("sample", drop=False)

samples_wc = [u.sample for u in samples.itertuples()]
subsamples_wc = [u.subsample for u in samples.itertuples()]

# All rule --> First one in a Snakemake pipeline.
rule all:
    input:
        filtered_graph=f"{OUTDIR}/PPI_network_filtered.tsv",
        annotation_file=f"{OUTDIR}/node_annotation_file.tsv"


# Define functions for rules.

def degs_input(wc):
    # Obtain input file from samples table.
    file = samples.loc[{wc.sample}]["file"]
    return file

def get_species(wc):
    # Obtain species information from samples table.
    file = samples.loc[{wc.sample}]["species"]
    return file


# Rule 1 reads the degs file for each dataset and sorts each by FDR and log2FC for downstream analysis.
# A filtered dataframe is given as an output.
rule filter_degs:
    input:
        degs_table=degs_input
    output:
        filtered_degs_table=f"{OUTDIR}/{{sample}}_filtered_degs_table.tsv"
    log: f"{LOGDIR}/{{sample}}_filter_degs.log"
    params:
        filter_table = config["parameters"]["filter_degs"]["filter_table"],
        outdir=OUTDIR,
        sample="{sample}",
        species=get_species
    conda: "envs/R.yaml"
    script:
        "scripts/step1_filter_degs.R"

# Rule 2 map the ligand receptor information to the DEGs table.
# A list of target genes (ligand & receptor) is produced.
# Also, the annotated DEGs table is exported.
rule map_ligand_receptor:
    input:
        filtered_degs_table=f"{OUTDIR}/{{sample}}_filtered_degs_table.tsv"
    output:
        LR_list=f"{OUTDIR}/{{sample}}_LR_list.txt",
        filtered_degs_table_ann=f"{OUTDIR}/{{sample}}_filtered_degs_table_LR_annotated.tsv"
    log: f"{LOGDIR}/{{sample}}_map_ligand_receptor.log"
    params:
        ligand_receptor_source = config["parameters"]["map_ligand_receptor"]["ligand_receptor_source"],
        outdir=OUTDIR,
        sample="{sample}"
    conda: "envs/R.yaml"
    script:
        "scripts/step2_map_ligand_receptor.R"

# Rule 3 filter the PPI database is filtered out according to the input genes.
# The network is exported as a TSV.
rule stringDB_PPI:
    input:
        LR_lists=expand("{OUTDIR}/{sample}_LR_list.txt", OUTDIR=OUTDIR, sample=samples_wc)
    output:
        graph=f"{OUTDIR}/PPI_network.tsv"
    log: f"{LOGDIR}/stringDB_PPI.log"
    params:
        outdir=OUTDIR,
        mode=config["parameters"]["stringDB_PPI"]["mode"],
        stringDB=config["parameters"]["stringDB_PPI"]["stringDB"]
    conda: "envs/R.yaml"
    script:
        "scripts/step3_stringDB_PPI.R"

# Rule 4 filter the PPI network according to tissue origin and ligand & receptor information.
# This fitlered PPI network is exported.
# Also the node (protein) information is saved.
rule filter_PPI:
    input:
        graph=f"{OUTDIR}/PPI_network.tsv",
        filtered_degs_tables_ann=expand("{OUTDIR}/{sample}_filtered_degs_table_LR_annotated.tsv", OUTDIR=OUTDIR, sample=samples_wc)
    output:
        filtered_graph=f"{OUTDIR}/PPI_network_filtered.tsv",
        annotation_file=f"{OUTDIR}/node_annotation_file.tsv",
    log: f"{LOGDIR}/filter_PPI.log"
    params:
        outdir=OUTDIR,
        samples=expand("{sample}", sample = samples_wc),
        PPI_score_threshold=config["parameters"]["filter_PPI"]["PPI_score_threshold"]
    conda: "envs/R.yaml"
    script:
        "scripts/step4_filter_PPI.R"
