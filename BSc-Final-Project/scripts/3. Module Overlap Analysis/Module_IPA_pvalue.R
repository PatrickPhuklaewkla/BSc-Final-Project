setwd("../Bsc-Final-Project/") # adjust to root directory
library(tidyverse)
library(data.table)
library(biomaRt)
library(rstatix)

source_dir  <- "data/modules/source"
output_dir  <- "data/modules/TST_module_pvalues/IPA"

module_files <- c(
  "ModuleA_URA_output.txt",
  "ModuleB_URA_output.txt",
  "ModuleC_URA_output.txt",
  "ModuleD_URA_output.txt",
  "IPA_DEG_A.txt",
  "IPA_DEG_B.txt",
  "IPA_DEG_C.txt",
  "IPA_DEG_D.txt"
)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir)

## Step 1: Load reference genes ####
listEnsembl()

ensembl <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

biomart <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name","gene_biotype"),
  mart = ensembl
)

type_protein <- c("protein_coding",
                  "IG_V_gene","IG_C_gene","IG_J_gene","IG_D_gene",
                  "TR_V_gene","TR_C_gene","TR_J_gene","TR_D_gene")

biomart_protein <- biomart %>%
  filter(gene_biotype %in% type_protein)

biomart_go <- getBM(
  attributes = c("hgnc_symbol"),
  filters = "go",
  values = c("GO:0005125"),
  mart = ensembl
)

cytokine_genes <- biomart_go %>%
  pull(hgnc_symbol) %>%
  unique()

## Step 2: Process module files iteratively ####

for (module_file in module_files) {
  
  input_path <- file.path(source_dir, module_file)
  
  df <- read.table(input_path,
                   header = TRUE,
                   sep = "\t",
                   skip = 2) %>%
    filter(B.H.corrected.p.value < 0.05) %>%
    filter(Upstream.Regulator %in% cytokine_genes)
  
  output_file <- paste0(
    tools::file_path_sans_ext(module_file),
    "_IPA_pvalue_output.csv"
  )
  
  write_csv(
    df,
    file.path(output_dir, output_file)
  )
}