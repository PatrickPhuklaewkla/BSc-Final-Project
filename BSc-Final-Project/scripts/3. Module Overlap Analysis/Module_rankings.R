setwd("../Bsc-Final-Project/") # adjust to root directory

library(tidyverse)
library(data.table)
library(biomaRt)
library(rstatix)

input_dir  <- "data/modules/TST_module_pvalues/collapsed"
output_dir <- "data/modules/TST_module_pvalues/ranking"
if (!dir.exists(output_dir)) dir.create(output_dir)

files <- list.files(input_dir, full.names = TRUE)

## ---- Filtering for genes of interest ----
genes_exact <- c("IFNG", "TNF", "TYPE_I_INTERFERON_COLLAPSED", "IFN_ALPHA_COLLAPSED")

for (f in files) {
  
  message("Processing: ", basename(f))
  
  df <- read_csv(f, show_col_types = FALSE) %>%
    select(regulator, p_adj, overlap_genes) %>% 
    filter(
      regulator %in% genes_exact |
        str_detect(regulator, "^IL17")
    )
  
  write_csv(df, file.path(output_dir, basename(f)))
}