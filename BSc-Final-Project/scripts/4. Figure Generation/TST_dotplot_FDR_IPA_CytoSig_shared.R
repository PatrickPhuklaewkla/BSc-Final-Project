setwd("../Bsc-Final-Project/data/")

library(tidyverse)

## Step 0: config settings ####
top_n <- NULL                     # NULL = all URAs
gradient_mode <- "global"       # "global" or "local"
rank_metric <- "adj_pvalue"     # "FDR" or "adj_pvalue"

ipa_file    <- "biomart_cytokine_restricted_IPA_URA_TSTd2vsSaline_fdrsig_648size.csv"
string_file <- "TST CytoSig/cytokine_restricted_cytosig_formatted_significant_UR_output.csv"

out_file <- paste0(
  "IPA_CytoSig_shared_dotplotURA_table_",
  rank_metric, "_",
  ifelse(is.null(top_n), "all", paste0("top", top_n)), "_",
  gradient_mode,
  ".pdf"
)

# Colours
high_col <- "#b30000"
low_col  <- "#1d4877"
missing_col <- "#F2F0EF"

## Step 1: Load files ####
ipa <- read_csv(ipa_file, show_col_types = FALSE)
string <- read_csv(string_file, show_col_types = FALSE)

## Step 2: Validate rank column ####
check_metric <- function(df, metric, label) {
  if (!metric %in% colnames(df)) {
    stop(paste(label, "missing column:", metric))
  }
}
check_metric(ipa, rank_metric, "IPA")
check_metric(string, rank_metric, "STRING")

## Step 3: Rank full tables ####
ipa_full <- ipa %>%
  arrange(.data[[rank_metric]], desc(abs(zscore))) %>%
  mutate(rank = row_number())

string_full <- string %>%
  arrange(.data[[rank_metric]], desc(abs(zscore))) %>%
  mutate(rank = row_number())

# optional top-n filter
ipa_use <- ipa_full
string_use <- string_full

if (!is.null(top_n)) {
  ipa_use <- slice_head(ipa_use, n = top_n)
  string_use <- slice_head(string_use, n = top_n)
}


max_len <- max(nrow(ipa_use), nrow(string_use))

## Step 4: Initialise palette ####
valid_ipa <- if (gradient_mode == "global") {
  ipa_full$regulator
} else if (gradient_mode == "local") {
  ipa_use$regulator
} else {
  stop("gradient_mode must be 'global' or 'local'")
}

valid_ipa <- as.character(valid_ipa)

ipa_palette <- colorRampPalette(c(high_col, "#f68838", "#fbb021", "#1b8a5a", low_col))(length(valid_ipa))
names(ipa_palette) <- valid_ipa

## Step 5: Build comparison table ####

ipa_ranks <- ipa_full %>%
  select(regulator, ipa_rank = rank)

string_ranks <- string_full %>%
  select(regulator, string_rank = rank)

paired_df <- inner_join(ipa_ranks, string_ranks, by = "regulator")
shared_regs <- paired_df$regulator

wilcox_res <- wilcox.test(
  paired_df$ipa_rank,
  paired_df$string_rank,
  paired = TRUE,
  alternative = "two.sided"
)

print(wilcox_res)

# ---- OPTIONAL TOP-N FILTER ----
if (!is.null(top_n)) {
  paired_df <- paired_df %>%
    filter(
      ipa_rank <= top_n | string_rank <= top_n
    )
}

# ---- REMOVE GAPS BY RE-RANKING ----
paired_df <- paired_df %>%
  arrange(ipa_rank) %>%
  mutate(ipa_rank = row_number()) %>%
  arrange(string_rank) %>%
  mutate(string_rank = row_number())

# ---- COLOUR LOGIC ----

paired_df <- paired_df %>%
  mutate(colour = ipa_palette[regulator])

# ---- LONG FORMAT FOR PAIRED PLOT ----

plot_df <- paired_df %>%
  pivot_longer(
    cols = c(ipa_rank, string_rank),
    names_to = "method",
    values_to = "rank"
  ) %>%
  mutate(
    x = ifelse(method == "ipa_rank", 1, 2)
  )

# ---- LABEL DATA ----

ipa_label_df <- plot_df %>%
  filter(
    method == "ipa_rank",
    regulator %in% shared_regs
  )

string_label_df <- plot_df %>%
  filter(
    method == "string_rank",
    regulator %in% shared_regs
  )

# ---- PAIRED DOT PLOT ----

p <- ggplot(paired_df) +
  
  geom_segment(
    aes(
      x = 1,
      xend = 2,
      y = ipa_rank,
      yend = string_rank,
      colour = colour
    ),
    linewidth = 0.8,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  
  geom_point(
    data = plot_df,
    aes(
      x = x,
      y = rank,
      fill = colour
    ),
    colour = "black",
    size = 3,
    shape = 21,
    na.rm = TRUE
  ) +
  
  # IPA labels (aligned column left of IPA nodes)
  geom_text(
    data = ipa_label_df,
    aes(
      x = 0.97,
      y = rank,
      label = regulator
    ),
    hjust = 1,
    size = 3,
    na.rm = TRUE
  ) +
  
  # STRING labels (aligned column right of STRING nodes)
  geom_text(
    data = string_label_df,
    aes(
      x = 2.03,
      y = rank,
      label = regulator
    ),
    hjust = 0,
    size = 3,
    na.rm = TRUE
  ) +
  
  scale_fill_identity() +
  scale_colour_identity() +
  
  scale_y_reverse(
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  
  scale_x_continuous(
    breaks = c(1,2),
    labels = c("IPA","CytoSig"),
    limits = c(0.6, 2.4)
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    y = paste0("Rank (1 = best ", rank_metric, ")"),
    title = paste0("Adjusted P-value Rank Comparison: Shared Regulators")
  )

# ---- SAVE ----

print(p)

ggsave(
  out_file,
  p,
  width = 6,
  height = max_len * 0.06,
  device = cairo_pdf
)
