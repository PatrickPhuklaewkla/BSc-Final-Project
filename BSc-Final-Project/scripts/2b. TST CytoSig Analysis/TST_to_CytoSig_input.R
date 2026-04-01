setwd("../data/")

library(biomaRt)
library(dplyr)
library(readr)

# read input file
df <- read_tsv("IPA_input_D2-TSTvsSaline.txt")

# remove Ensembl version numbers if present
df$Id <- sub("\\..*", "", df$Id)

# connect to Ensembl
mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

# map Ensembl -> gene symbol
mapping <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = df$Id,
  mart = mart
)

# merge mapping
df2 <- df %>%
  left_join(mapping, by = c("Id" = "ensembl_gene_id")) %>%
  filter(hgnc_symbol != "")

# compute CytoSig score (log2FoldChange)
df2 <- df2 %>%
  mutate(score = log2FoldChange)

# handle duplicated gene symbols
df2 <- df2 %>%
  group_by(hgnc_symbol) %>%
  summarise(score = max(score), .groups = "drop")

# convert to named vector for CytoSig
cytosig_input <- df2$score
names(cytosig_input) <- df2$hgnc_symbol

# save
write.table(
  df2[,c("hgnc_symbol","score")],
  "cytosig_input.txt",
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)
