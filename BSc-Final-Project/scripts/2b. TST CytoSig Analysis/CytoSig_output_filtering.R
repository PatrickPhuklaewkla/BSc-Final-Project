setwd("../data/")

library(tidyverse)
library(data.table)
library(biomaRt)
library(rstatix)

## Load cytokine activity terms for filtering ####
listEnsembl() # check current version
#listEnsemblArchives() # check archived versions

ensembl <- useEnsembl(biomart = "genes",
                      #version = 114, # only needed if using archived version
                      dataset = "hsapiens_gene_ensembl")

attributes <- listAttributes(ensembl)

biomart <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype"),
  mart = ensembl
)

# reference = ~20k protein-coding genes
type <- unique(biomart$gene_biotype)
type_protein <- c("protein_coding",
                  "IG_V_gene","IG_C_gene","IG_J_gene","IG_D_gene",
                  "TR_V_gene","TR_C_gene","TR_J_gene","TR_D_gene")
biomart_protein <- biomart %>%
  filter(gene_biotype %in% type_protein)

# Cytokine (/ cytokine receptor genes) GO query
biomart_go <- getBM(
  attributes = c("hgnc_symbol"),
  filters = "go",
  values = c("GO:0005125"), # "GO:0005125": Cytokine Activity, "GO:0004896": Cytokine Receptor Activity
  mart = ensembl
)

cytokine_genes <- biomart_go %>%
  dplyr::pull(hgnc_symbol) %>%
  unique()

## Load CIE output file ####
cytosig <- read.csv("data/cytosig_output.csv")

## Filter CIE file ####
cytosig <- cytosig %>%
  filter(Pvalue < 0.05) %>%
  filter(Zscore > 0) %>%
  filter(ID %in% cytokine_genes)

## Write to file ####
write_csv(cytosig,
          "data/cytokine_restricted_cytosig_significant_UR_output.csv")
