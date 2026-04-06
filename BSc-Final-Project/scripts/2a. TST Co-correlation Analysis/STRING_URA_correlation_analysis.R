## Identification of statistically significant co-regulated upstream regulator modules (TST_D2 vs saline)

# script adapted from Blanca Sanz-Magallon Duque De Estrada

setwd("../data/")

library(tidyverse)

## A. Input files - replace as necessary
tpm <- read.csv("tpm_D2-transcriptome.csv", row.names = 1)
ipa <- read.csv("cytokine_restricted_STRING_significant_UR_output.csv", header = TRUE) # map STRING output to established ipa analysis
meta <- read.csv("RNAseq_metadata.csv", header = TRUE, row.names = 1)

## B. Create a co-correlation matrix of the TST transcriptome (including only relevant TST samples)
meta.ss <- meta %>% filter(Stimulant %in% c("TST_D2"))
tpm.ss <- tpm[,rownames(meta.ss)]

TST_transcriptome <- tpm.ss
TST_transcriptome <- TST_transcriptome[order(rownames(TST_transcriptome)),]
r_tst <- cor(t(TST_transcriptome), method = "spearman")

## C. Format the IPA output file
colnames(ipa) <- c("regulator","FDR","count","target")
ipa <- tibble::rowid_to_column(ipa, "ID") #add an ID column corresponding to each cluster
#ipa$target <- as.character(ipa$target) #target column is integer -> convert to character for next step
ipa_sep <- separate_rows(ipa, target, sep = ",|/", convert = TRUE) #separate clusters into individual rows
ipa_sep_reduced <- subset(ipa_sep, target %in% colnames(r_tst)) #only the genes found in the TST transcriptome co-correlation matrix
ipa_sep_reduced <- ipa_sep_reduced %>%
  relocate(count, .after = last_col())

## D. Calculate the average of all pairwise correlations for the target genes of each upstream regulator cluster
avg_interactomes <- matrix(nrow=0, ncol=2) # empty data frame to add the calculated average correlation values
colnames(avg_interactomes) <- c("regulator", "avg_corr")
for (i in 1:max(ipa_sep_reduced$ID)){
  hub <- ipa_sep_reduced %>% filter(ID == i)
  genes <- hub$target[which(hub$target %in% row.names(r_tst))]
  r_tst_network <- r_tst[genes,genes]
  lower_matrix <- lower.tri(r_tst_network)
  r_tst_network_lower <- r_tst_network[lower_matrix]
  mean_interactome <- mean(r_tst_network_lower[r_tst_network_lower>=0])
  avg_interactomes <- rbind(avg_interactomes, c(as.character(unique(hub$regulator)), mean_interactome))
}
# add the calculated cluster average correlations to a data frame containing the other information
ipa_sep_reduced$target <- NULL
ipa_sep_reduced <- unique(ipa_sep_reduced)
ipa_sep_reduced <- merge(ipa_sep_reduced, avg_interactomes, by.x = "regulator", by.y = "regulator") # merge the ipa output file and the correlation values by common regulator name

## E. Generate frequency distributions of the expression of randomly selected groups of genes of a range of sizes
## notes:
##   p = the maximum size of the group of genes to test (ie. range 4 - k). 4 is the lowest number for the range, the correlation function doesn't work for size < 4
##   r = the number of random samples used to calculate the frequency distribution for a given size -- the larger the number the more accurate, but at the cost of longer running time
p <- max(ipa_sep_reduced$count) # check ipa_sep_reduced for largest count
r <- 100
gene_names <- rownames(TST_transcriptome)

random_distributions_100 <- matrix(nrow = 0, ncol = 6) # empty matrix to fill with the distribution values
colnames(random_distributions_100) <- c("size", "mean", "sd", "84.13%", "97.72%",	"99.87%")
for (k in (4:p)){
  random_genes <- matrix(nrow = 0, ncol = 2)
  colnames(random_genes) <- c("iteration", "average_correlation")
  for (i in (1:r)){ # Do this 100 times (could also try doing 1000 times if it doesn't take too long for the loop to run)
    genes <- sample(gene_names, k) # random sample from TST transcriptome (note: no seed set, so repeat runs will yield slightly different results)
    r_tst_network <- r_tst[genes,genes]
    lower_matrix <- lower.tri(r_tst_network)
    r_tst_network_lower <- r_tst_network[lower_matrix]
    mean_interactome <- mean(r_tst_network_lower[r_tst_network_lower>=0]) # Calculate the average correlation for each of the 100 (or 1000) random clusters
    random_genes <- rbind(random_genes, c(i, mean_interactome))
  }
  mean <- mean(as.numeric(random_genes[,2]))
  sd <- sd(as.numeric(random_genes[,2]))
  random_distributions_100 <- rbind(random_distributions_100, c(k, mean, sd, mean + sd, mean + 2*sd, mean + 3*sd))
}
colnames(random_distributions_100) <- c("size", "mean", "sd", "X84.13", "X97.72",	"X99.86")
random_distributions_100_df <- as.data.frame(random_distributions_100)

## F. Determining which clusters are FDR significant
ipa_avg_corr <- ipa_sep_reduced[ipa_sep_reduced$count >= 4,] # take out clusters smaller than 4
ipa_avg_corr$rc_mean <- random_distributions_100_df$mean[match(ipa_avg_corr$count, random_distributions_100_df$size)] #mean
ipa_avg_corr$rc_sd <- random_distributions_100_df$sd[match(ipa_avg_corr$count, random_distributions_100_df$size)] #sd
ipa_avg_corr$avg_corr <- as.numeric(ipa_avg_corr$avg_corr)
ipa_avg_corr$zscore <- (ipa_avg_corr$avg_corr - ipa_avg_corr$rc_mean)/ipa_avg_corr$rc_sd #zscores
ipa_avg_corr$pvalue <- pnorm(ipa_avg_corr$zscore, lower.tail = FALSE) #p-values
ipa_less_0.05 <- ipa_avg_corr[which(ipa_avg_corr$pvalue <= 0.05),] # p-values <= 0.05
fdr_threshold <- nrow(ipa_less_0.05)
ipa_less_0.05$adj_pvalue <- ipa_less_0.05$pvalue*fdr_threshold # calculating the adjusted p-value
ipa_significant <- ipa_less_0.05[which(ipa_less_0.05$adj_pvalue <= 0.05),] # selecting clusters with adj. p-value <= 0.05

## Add type column
ipa_significant$type <- "cytokine"

ipa_significant <- ipa_significant %>%
  dplyr::relocate(type, .after = 2)

## G. Output file
write.csv(ipa_significant, paste0("cytokine_restricted_STRING_URA_TSTd2vsSaline_fdrsig_",p,"size.csv"), row.names = F) # clusters with adj. p-value <= 0.05
write.csv(random_distributions_100, paste0("random_distrib_4-",p,"size_100its_tpm_TST_avg-correl.csv"))  # change file name as appropriate (inc. iterations and range)
write.csv(ipa_avg_corr, "ipa_avg_corr.csv", row.names = FALSE)

#random_distributions_100 <- read.csv("random_distrib_4-600size_100its_tpm_TST_avg-correl.csv")

## H. Plotting average cluster expression vs. frequency distribution
library(ggplot2)

ipa_avg_corr <- read.csv("ipa_avg_corr.csv")
random_distributions_100 <- read.csv("random_distrib_4-1309size_100its_tpm_TST_avg-correl.csv") # read in output file from Script 2 if not running steps A-D above
random_distributions_100_df <- as.data.frame(random_distributions_100) #plotting only works with data frame
ipa_significant <- read.csv("cytokine_restricted_STRING_URA_TSTd2vsSaline_fdrsig_1309size.csv", header = TRUE)

# load files if not running steps A-C above
#ipa_avg_corr <- read.csv("tpm_TST_avg-correl.csv") # output file from Script 1


# Cluster correlations with frequency distributions plot
random_genes_plot <- ggplot(ipa_avg_corr, aes(avg_corr, count)) + 
  xlab("Average correlation coefficient") + ylab("Network size (number of target genes)") + 
  #ggtitle("Average Cluster Correlation Coefficient") +
  geom_point(data = random_distributions_100_df, aes(X97.72, size), color = "palegreen3") + 
  geom_point(data = random_distributions_100_df, aes(X84.13, size), color = "steelblue2") +
  geom_point(data = random_distributions_100_df, aes(X99.86, size), color = "sienna2") +
  geom_point(color = "darkgrey")+
  theme_light(base_size = 14)
plot(random_genes_plot)

ipa_avg_corr$significant <- ipa_avg_corr$regulator %in% ipa_significant$regulator
ggplot(ipa_avg_corr, aes(avg_corr, count, color = significant)) +
  xlab("Average correlation coefficient") + ylab("Network size (number of target genes)") +
  ggtitle("Average Cluster Correlation Coefficient") +
  geom_point() +
  geom_point(data=random_distributions_100_df, aes(X97.72, size), color="palegreen3") +
  geom_point(data=random_distributions_100_df, aes(X84.13, size), color="steelblue2") +
  geom_point(data=random_distributions_100_df, aes(X99.86, size), color="sienna2") +
  scale_color_manual(values=c("grey60","blue")) +
  theme_light(base_size = 14)

# Plotting the distributions only
random_distributions_plot <- ggplot(random_distributions_100_df, aes(X97.72, size)) + 
  geom_point(color = "palegreen3") +
  xlab("Average correlation coefficient") + ylab("Network size (number of target genes)") + #replace axes names here
  #ggtitle("Average Cluster Correlation Coefficient") + # replace title here
  geom_point(data = random_distributions_100_df, aes(X84.13, size), color = "steelblue2") +
  geom_point(data = random_distributions_100_df, aes(X99.86, size), color = "sienna2")+
  theme_light(base_size = 14)
plot(random_distributions_plot)
