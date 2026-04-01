setwd("../Bsc-Final-Project/") # adjust to root directory

library(tidyverse)
library(data.table)
library(biomaRt)

## Step 0: Directories ####

source_dir  <- "data/modules/source"
output_dir  <- "data/modules/TST_module_pvalues/CytoSig"
cytosig_file <- "data/cytosig/signature.centroid.expand"

module_files <- c(
  "Module A.csv",
  "Module B.csv",
  "Module C.csv",
  "Module D.csv",
  "DEG_A_clean.csv",
  "DEG_B_clean.csv",
  "DEG_C_clean.csv",
  "DEG_D_clean.csv"
)

## Step 1: Load reference gene universe (Ensembl protein coding) ####

ensembl <- useEnsembl(
  biomart = "genes",
  #version = 115,
  dataset = "hsapiens_gene_ensembl"
)

biomart <- getBM(
  attributes = c("ensembl_gene_id","external_gene_name","gene_biotype"),
  mart = ensembl
)

type_protein <- c(
  "protein_coding",
  "IG_V_gene","IG_C_gene","IG_J_gene","IG_D_gene",
  "TR_V_gene","TR_C_gene","TR_J_gene","TR_D_gene"
)

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


ref <- biomart_protein %>%
  pull(external_gene_name) %>%
  unique()

## Step 2: Load CytoSig cytokine signature matrix ####

# Matrix structure:
# rows = genes
# columns = cytokines
# values = cytokine response coefficients

# Read header line
header <- readLines(cytosig_file, n = 1)
header <- strsplit(header, "\t")[[1]]

# Add missing first column name
header <- c("gene", header)
cytosig <- read.delim(
  cytosig_file,
  sep = "\t",
  header = FALSE,
  skip = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Rename first column as gene symbol column
colnames(cytosig) <- header

## Convert cytokine matrix into long format ####

cytosig_long <- cytosig %>%
  pivot_longer(
    cols = -gene,
    names_to = "regulator",
    values_to = "score"
  )

# Define cytokine target genes
# (top 200 strongest responding genes per cytokine)

cytosig_targets <- cytosig_long %>%
  group_by(regulator) %>%
  slice_max(order_by = abs(score), n = 200) %>%
  ungroup()


## Convert to regulator → gene list format ####

cytosig_list <- cytosig_targets %>%
  group_by(regulator) %>%
  summarise(target = list(unique(gene))) %>%
  deframe()

cytosig_list <- lapply(cytosig_list, function(x){
  data.frame(target = x)
})

## Step 4: Overlap function (Fisher enrichment test) ####
compute_overlap <- function(ur_df, tst_genes, ref){
  
  ur_targets <- unique(ur_df$target)
  tst <- unique(tst_genes)
  
  overlap_genes <- intersect(tst, ur_targets)
  
  a <- length(overlap_genes)
  b <- length(setdiff(tst, ur_targets))
  c <- length(setdiff(ur_targets, tst))
  d <- length(setdiff(ref, union(tst, ur_targets)))
  
  tbl <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
  
  list(
    table = tbl,
    overlap_genes = overlap_genes
  )
}

## Step 5: Loop through modules ####

all_results <- list()

for(file in module_files){
  
  module_name <- tools::file_path_sans_ext(file)
  
  message("Processing ", module_name)
  
  # Load module gene list (robust to header / no header)
  module_genes <- read.csv(
    file.path(source_dir, file),
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # ensure first column is treated as gene symbol column
  colnames(module_genes)[1] <- "external_gene_name"
  
  module_genes <- module_genes %>%
    mutate(external_gene_name = trimws(external_gene_name)) %>%
    separate_rows(external_gene_name, sep = "///") %>%
    mutate(external_gene_name = trimws(external_gene_name)) %>%
    filter(
      !is.na(external_gene_name),
      external_gene_name != ""
    ) %>%
    distinct() %>%
    filter(external_gene_name %in% ref)
  
  tst <- module_genes$external_gene_name
  
  # Compute overlap against all cytokine signatures
  
  overlap_tables <- lapply(
    cytosig_list,
    compute_overlap,
    tst_genes = tst,
    ref = ref
  )
  
  overlap_tables_f <- overlap_tables[
    sapply(overlap_tables, function(x) length(x$overlap_genes)) > 0
  ]
  
  fisher_results <- lapply(overlap_tables_f, function(x){
    fisher.test(x$table, alternative="greater")
  })
  
  # Assemble results
  
  results_df <- tibble(
    regulator = names(overlap_tables_f),
    p_value = sapply(fisher_results, `[[`, "p.value"),
    overlap_genes = sapply(
      overlap_tables_f,
      function(x) paste(sort(x$overlap_genes), collapse=",")
    ),
    n_overlap = sapply(
      overlap_tables_f,
      function(x) length(x$overlap_genes)
    )
  )
  
  results_df_sig <- results_df %>%
    mutate(fdr = p.adjust(p_value, method="BH")) %>%
    filter(fdr < 0.05) %>%
    #filter(regulator %in% cytokine_genes) %>%
    arrange(fdr)
  
  write_csv(
    results_df_sig,
    file.path(
      output_dir,
      paste0(module_name, "_cytosig_enrichment.csv")
    )
  )
  
  all_results[[module_name]] <- results_df_sig
}
