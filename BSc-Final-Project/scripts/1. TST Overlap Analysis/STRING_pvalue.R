setwd("../Bsc-Final-Project/") # adjust to root directory

library(tidyverse)
library(data.table)
library(biomaRt)
library(rstatix)

## Step 1: Load reference genes ####
# https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/accessing_ensembl.html

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

## Step 2: Load interaction database and subset to reference genes ####
# STRING: downloaded from https://string-db.org/cgi/download?sessionId=bgFsQUBM6MYS&species_text=Homo+sapiens
string_db <- fread("interaction_db_files/9606.protein.links.full.v12.0.txt.gz") 
string_proteins <- fread("interaction_db_files/9606.protein.info.v12.0.txt.gz") # protein annotations, eg. common names

prot_short <- string_proteins %>%
  dplyr::rename(protein1 = "#string_protein_id") %>%
  mutate(protein2 = protein1,
         name1 = preferred_name,
         name2 = preferred_name) %>%
  dplyr::select(protein1,name1,protein2,name2)
prot1 <- prot_short %>% dplyr::select(protein1,name1)
prot2 <- prot_short %>% dplyr::select(protein2,name2)
string_db_annot <- string_db %>% left_join(prot1) %>% left_join(prot2) 

diff1 <- setdiff(string_db_annot$name1,biomart_protein$external_gene_name)
diff2 <- setdiff(string_db_annot$name2,biomart_protein$external_gene_name)
all(diff1 %in% diff2) # TRUE

string_db_final <- string_db_annot %>%
  filter(name1 %in% biomart_protein$external_gene_name &
           name2 %in% biomart_protein$external_gene_name)

## Step 3: Select relevant fields for all upstream regulators and tibble data frame ####
ur_network_string <- string_db_final %>%
  dplyr::select(name1,name2,combined_score) %>%
  unique()

ur_network_string_sort <- ur_network_string %>% arrange(name1)

ur_list <- ur_network_string_sort %>%
  group_by(name1) %>%
  group_split()

ur_names <- ur_network_string_sort %>%
  pull(name1) %>%
  unique()

names(ur_list) <- ur_names

# network_sizes <- sapply(ur_list, nrow) # optional for checking summary statistics
# summary(network_sizes)

## Step 4: Load gene list of interest (genes up-regulated in D2 TST) and subset to reference genes ####

tst <- read.table("data/IPA_input_D2-TSTvsSaline.txt",
                  header = TRUE)
tst_annot <- tst %>%
  dplyr::rename(ensembl_gene_id = "Id") %>%
  left_join(biomart_protein) %>%
  dplyr::select(external_gene_name) %>%
  filter(!is.na(external_gene_name) &
           !external_gene_name == "") %>%
  unique()

## Step 5: calculate overlap p-value for every upstream regulator ####

ref <- biomart_protein %>% pull(external_gene_name) %>% unique()

# Define function
compute_overlap <- function(ur_df, tst_annot, ref) {
  
  ur_targets <- unique(ur_df$name2)
  tst <- unique(tst_annot$external_gene_name)
  
  # a = overlapping genes used in Fisher test
  overlap_genes <- intersect(tst, ur_targets)
  
  a <- length(overlap_genes)
  b <- length(setdiff(tst, ur_targets))
  c <- length(setdiff(ur_targets, tst))
  d <- length(setdiff(ref, union(tst, ur_targets)))
  
  tbl <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  list(
    table = tbl,
    overlap_genes = overlap_genes
  )
}
# Calculate overlap p-values
overlap_tables <- lapply(
  ur_list,
  compute_overlap,
  tst_annot = tst_annot,
  ref = ref
)

fisher_results <- lapply(overlap_tables, function(x) {
  fisher.test(x$table, alternative = "greater") 
})

## Step 6: Fetch each UR and associated fields, correcting for multiple testing using Benjamini-Hochberg Procedure ####

results_df <- tibble(
  regulator = names(overlap_tables),
  p_value = sapply(fisher_results, `[[`, "p.value"),
  overlap_genes = sapply(overlap_tables, function(x) # collapses target list into an IPA-like format
    paste(sort(x$overlap_genes), collapse = ",")
  ),
  n_overlap = sapply(overlap_tables, function(x)
    length(x$overlap_genes)
  )
)

results_df_sig <- results_df %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(n_overlap > 0, p_adj < 0.05) %>%
  filter(regulator %in% cytokine_genes) %>%   # cytokine(/cytokine receptor) only
  dplyr::select(regulator, p_adj, n_overlap, overlap_genes) %>%
  arrange(p_adj)

write_csv(results_df_sig,
          "data/cytokine_restricted_STRING_significant_UR_output.csv")
