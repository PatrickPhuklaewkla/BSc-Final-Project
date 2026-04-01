setwd("../Bsc-Final-Project/") # adjust to root directory

library(tidyverse)
library(data.table)
library(biomaRt)
library(rstatix)

## Step 0: Initialise interferon definitions ####
type1_ifn_genes <- c(
  "IFNA1","IFNA2","IFNA4","IFNA5","IFNA6","IFNA7","IFNA8",
  "IFNA10","IFNA13","IFNA14","IFNA16","IFNA17","IFNA21",
  "IFNB1","IFNE","IFNK","IFNW1"
)

ifn_alpha_genes <- grep("^IFNA", type1_ifn_genes, value = TRUE)

collapse_family <- function(df, genes, label) {
  fam_rows <- df %>% filter(regulator %in% genes)
  if (nrow(fam_rows) == 0) return(NULL)
  
  combined_genes <- fam_rows$overlap_genes %>%
    strsplit(",") %>%
    unlist() %>%
    unique() %>%
    sort()
  
  tibble(
    regulator = label,
    p_adj = min(fam_rows$p_adj, na.rm = TRUE),
    n_overlap = length(combined_genes),
    overlap_genes = paste(combined_genes, collapse = ",")
  )
}

## Step 1: Iterate over STRING + IPA files ####

input_dirs <- c(
  "data/modules/TST_module_pvalues/STRING",
  "data/modules/TST_module_pvalues/IPA"
)

output_dir <- "data/modules/TST_module_pvalues/collapsed"
if (!dir.exists(output_dir)) dir.create(output_dir)

# collect files from both directories
all_files <- unlist(lapply(input_dirs, list.files, full.names = TRUE))

# keep STRING or IPA files
target_files <- all_files[grepl("STRING|IPA", basename(all_files), ignore.case = TRUE)]

for (f in target_files) {
  
  message("Processing: ", basename(f))
  
  df <- read_csv(f, show_col_types = FALSE)
  
  # IPA catch-case: harmonise column names
  if (grepl("IPA", basename(f), ignore.case = TRUE)) {
    
    rename_map <- c(
      "Upstream.Regulator" = "regulator",
      "B.H.corrected.p.value" = "p_adj",
      "Target.Molecules.in.Dataset" = "overlap_genes"
    )
    
    existing <- intersect(names(rename_map), colnames(df))
    if (length(existing) > 0) {
      df <- df %>% rename(!!!setNames(existing, rename_map[existing]))
    }
  }
  
  if (all(c("regulator", "p_adj", "overlap_genes") %in% colnames(df))) {
    
    collapsed_ifn_summary <- bind_rows(
      collapse_family(df, type1_ifn_genes, "TYPE_I_INTERFERON_COLLAPSED"),
      collapse_family(df, ifn_alpha_genes, "IFN_ALPHA_COLLAPSED")
    )
    
    df <- df %>%
      bind_rows(collapsed_ifn_summary) %>%
      arrange(p_adj)
  } else {
    warning("Skipping IFN collapsing (missing required columns): ", basename(f))
  }
  
  write_csv(df, file.path(output_dir, basename(f)))
}