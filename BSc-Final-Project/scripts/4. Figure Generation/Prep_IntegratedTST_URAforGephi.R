
# Prep packages & themes --------------------------------------------------

library(data.table)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(pROC)
library(gtsummary)
library(magrittr)

#My_Theme
t = 10 #size of text
m = 5 #size of margin around text
tc = "black" #colour of text
My_Theme = theme(
  axis.title.x = element_text(size = t, face = "bold", margin = margin(t = m)),
  axis.text.x = element_text(size = t, face = "bold", colour = tc, angle = 0, hjust = 0.5),
  axis.title.y.left = element_text(size = t, face = "bold", margin = margin(r = m)),
  axis.title.y.right = element_text(size = t, face = "bold", margin = margin(l = m)),
  axis.text.y = element_text(size = t, face = "bold", colour = tc),
  legend.title = element_text(size=t, face = "bold", colour = tc),
  legend.text = element_text(size=t, face = "bold", colour = tc),
  plot.title = element_text(size=t, face = "bold", colour = tc),
  strip.text = element_text(size=t, face = "bold", colour = tc),
  strip.background = element_rect(fill = "gray90", colour = "black", size = 0.5),
  panel.border = element_rect(fill = NA, linewidth = 0.5, colour = tc),
  panel.background = element_rect(fill = "gray97"),
  plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  legend.position = "right", legend.justification = "top"
)

cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73",
               "#F0E442","#0072B2","#D55E00","#CC79A7")

# Select files ------------------------------------------------------------

setwd("C:/Users/patri/OneDrive - University College London/Documents/Research Project Code/Signor/data")

URA_raw <- fread("cytokine_restricted_STRING_URA_TSTd2vsSaline_fdrsig_1309size.csv")

UR_targets <- fread("cytokine_restricted_STRING_significant_UR_output.csv") %>%
  select(regulator, overlap_genes)

URA <- URA_raw %>%
  select(regulator, zscore) %>%                # keep needed columns
  left_join(UR_targets, by = "regulator") %>%  # attach new targets
  mutate(type = "cytokine") %>%                # define regulator class
  select(regulator, overlap_genes, type, zscore) %>%
  set_colnames(c(
    "Source",
    "Target",
    "Class",
    "Zscore"
  ))


edges <- URA %>%
  distinct(Target, .keep_all = T) %>% 
  separate_rows(Target, sep = ",|/", convert = TRUE) %>% 
  mutate(Type="Directed") %>% 
  relocate(Type, .before=Class)

regulator <- edges %>% 
  select(Source,Class,Zscore) %>%
  rename(Id=Source)

target <- edges%>% 
  select(Target) %>% 
  rename(Id=Target) %>% 
  mutate(Class="target") %>% 
  mutate(Zscore=1)

nodes <- rbind(regulator,target) %>% 
  mutate(Label=Id) %>% 
  select(Id,Label,Class,Zscore) %>%
  distinct(Id,.keep_all = T)

write_tsv(nodes, "IntegratedTST_URA_GephiNodes.tsv")
write_tsv(edges, "IntegratedTST_URA_GephiEdges.tsv")

