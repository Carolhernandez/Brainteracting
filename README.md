# Brainteracting
Pipeline for dissecting metastasis-microenvironment crosstalk by gene expression interacting pairs.

**Setup**
The setup of the pipeline consists in the modification of three configuration files, indicating the desired parameters and the location of the input files.
A general description of these files follows. See the Usage section for more details.

**Configuration files**


config.yaml contains all pipeline parameters.

samples.tsv contains information on the samples to be analysed.

units.tsv: contains information on the different data files associated to the samples to be analysed.

**Input files**
Differential expression tables from both compartents:
DEGs_microenvironment.csv
DEGs_metastasis.csv
