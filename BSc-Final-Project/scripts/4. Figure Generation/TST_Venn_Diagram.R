setwd("../Bsc-Final-Project/data")

library(tidyverse)
library(VennDiagram)
library(grid)

# Load data
ipa <- read_csv("biomart_cytokine_restricted_IPA_URA_TSTd2vsSaline_fdrsig_648size.csv",
                show_col_types = FALSE)

string <- read_csv("cytokine_restricted_STRING_URA_TSTd2vsSaline_fdrsig_1309size.csv",
                   show_col_types = FALSE)

# Unique regulator sets
ipa_regs <- unique(ipa$regulator)
string_regs <- unique(string$regulator)

# Counts
area1 <- length(ipa_regs)
area2 <- length(string_regs)
cross_area <- length(intersect(ipa_regs, string_regs))

# Draw horizontal Venn
venn.plot <- draw.pairwise.venn(
  area1 = area1,
  area2 = area2,
  cross.area = cross_area,
  category = c("IPA", "STRING"),
  fill = c("#E69F00", "#009E73"),   # Orange, Green
  alpha = 0.6,
  cat.pos = c(-90, 90),             # Adjust label position
  cat.dist = c(0.05, 0.05),
  cex = 1.5,
  cat.cex = 1.4,
  scaled = TRUE,                     # Proportional areas
  margin = 0.06,
  ext.dist = 0.05,
  ext.length = 0.9
)

# Save as vector PDF
pdf("URA_IPA_vs_STRING_Venn_Diagram.pdf", width = 7, height = 5)
grid.draw(venn.plot)
dev.off()
