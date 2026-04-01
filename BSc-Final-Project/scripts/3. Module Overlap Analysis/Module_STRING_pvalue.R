setwd("../Bsc-Final-Project/") # adjust to root directory

library(tidyverse)
library(data.table)
library(biomaRt)
library(rstatix)

source_dir  <- "data/modules/source"
output_dir  <- "data/modules/TST_module_pvalues/STRING"

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

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(output_dir)) dir.create(output_dir)

## Step 1: Load reference genes ####
listEnsembl()

ensembl <- useEnsembl(
  biomart = "genes",
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

biomart_go <- getBM(
  attributes = c("hgnc_symbol"),
  filters = "go",
  values = c("GO:0005125"),
  mart = ensembl
)

cytokine_genes <- biomart_go %>%
  pull(hgnc_symbol) %>%
  unique()

## Step 2: Load STRING interaction database ####
string_db <- fread("interaction_db_files/9606.protein.links.full.v12.0.txt.gz")
string_proteins <- fread("interaction_db_files/9606.protein.info.v12.0.txt.gz")

prot_short <- string_proteins %>%
  rename(protein1 = "#string_protein_id") %>%
  mutate(
    protein2 = protein1,
    name1 = preferred_name,
    name2 = preferred_name
  ) %>%
  select(protein1,name1,protein2,name2)

prot1 <- prot_short %>% select(protein1,name1)
prot2 <- prot_short %>% select(protein2,name2)

string_db_annot <- string_db %>%
  left_join(prot1) %>%
  left_join(prot2)

string_db_final <- string_db_annot %>%
  filter(
    name1 %in% biomart_protein$external_gene_name &
      name2 %in% biomart_protein$external_gene_name
  )

## Step 3: Build UR network ####
ur_network_string <- string_db_final %>%
  select(name1,name2,combined_score) %>%
  unique()

ur_network_string_sort <- ur_network_string %>%
  arrange(name1)

ur_list <- ur_network_string_sort %>%
  group_by(name1) %>%
  group_split()

ur_names <- ur_network_string_sort %>%
  pull(name1) %>%
  unique()

names(ur_list) <- ur_names

ref <- biomart_protein %>%
  pull(external_gene_name) %>%
  unique()

## Step 4: Define overlap analysis function ####
compute_overlap <- function(ur_df, tst_annot, ref) {
  
  ur_targets <- unique(ur_df$name2)
  tst <- unique(tst_annot$external_gene_name)
  
  overlap_genes <- intersect(tst, ur_targets)
  
  a <- length(overlap_genes)
  b <- length(setdiff(tst, ur_targets))
  c <- length(setdiff(ur_targets, tst))
  d <- length(setdiff(ref, union(tst, ur_targets)))
  
  tbl <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)
  
  list(
    table = tbl,
    overlap_genes = overlap_genes
  )
}

## Step 5: Iterate overlap analysis over all module files ####
for (module_file in module_files) {
  
  message("Processing: ", module_file)
  
  input_path <- file.path(source_dir, module_file)
  
  module_genes <- read.csv(
    file.path(source_dir, module_file),
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  # ensure first column is gene symbol column
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
  
  cytokine_sig_annot <- module_genes
  tst <- module_genes$external_gene_name
  
  overlap_tables <- lapply(
    ur_list,
    compute_overlap,
    tst_annot = cytokine_sig_annot,
    ref = ref
  )
  
  overlap_tables_f <- overlap_tables[
    sapply(overlap_tables, function(x) length(x$overlap_genes)) > 0
  ]
  
## Step 6: Fetch each UR and associated p value, correcting for multiple testing using Benjamini-Hochberg Procedure and write to file ####
  
  fisher_results <- lapply(overlap_tables_f, function(x) {
    fisher.test(x$table, alternative = "greater")
  })
  
  results_df <- tibble(
    regulator = names(overlap_tables_f),
    p_value = sapply(fisher_results, `[[`, "p.value"),
    overlap_genes = sapply(
      overlap_tables_f,
      function(x) paste(sort(x$overlap_genes), collapse = ",")
    ),
    n_overlap = sapply(
      overlap_tables_f,
      function(x) length(x$overlap_genes)
    )
  )
  
  results_df_sig <- results_df %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    filter(p_adj < 0.05) %>%
    filter(regulator %in% cytokine_genes) %>%
    filter(n_overlap > 0) %>%
    select(regulator, p_adj, n_overlap, overlap_genes) %>%
    arrange(p_adj)
  
  output_file <- paste0(
    tools::file_path_sans_ext(module_file),
    "_STRING_pvalue_output.csv"
  )
  
  write_csv(
    results_df_sig,
    file.path(output_dir, output_file)
  )
}