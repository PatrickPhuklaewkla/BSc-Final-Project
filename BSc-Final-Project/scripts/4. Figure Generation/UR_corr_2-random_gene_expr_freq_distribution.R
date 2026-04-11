########## 2. Random gene frequency distributions and plots ########## Blanca, October 8th, 2020
library(dplyr)

## A. Input files - replace as necessary
tpm <- read.csv("../tpm_TST-transcriptome_Wilcox_GeneSymbol_2022.csv", na.strings = "")
tpm <- tpm %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames("GeneSymbol")
meta <- read.csv("../../Akshay/data/metadata.csv", header = TRUE, row.names = 2)
ipa_sep_reduced <- read.csv("tpm_TST-2022_avg-correl.csv") # Output file from Script 1

## B. Create a co-correlation matrix of the TST transcriptome (including only relevant TST samples)
meta.ss <- subset(meta, Stimulant == "TST_D7")
tpm.ss <- tpm[,rownames(meta.ss)]

TST_transcriptome <- tpm.ss
TST_transcriptome <- TST_transcriptome[order(rownames(TST_transcriptome)),]
r_tst <- cor(t(TST_transcriptome), method = "spearman")

## C. Generate frequency distributions of the expression of randomly selected groups of genes of a range of sizes
## notes:
##   p = the maximum size of the group of genes to test (ie. range 4 - k). 4 is the lowest number for the range, the correlation function doesn't work for size < 4
##   r = the number of random samples used to calculate the frequency distribution for a given size -- the larger the number the more accurate, but at the cost of longer running time
p <- 570 # check ipa_sep_reduced for largest count
r <- 100
gene_names <- rownames(TST_transcriptome)

random_distributions_100 <- matrix(nrow = 0, ncol = 6) # empty matrix to fill with the distribution values
colnames(random_distributions_100) <- c("size", "mean", "sd", "84.13%", "97.72%",	"99.87%")
for (k in (4:p)){
  random_genes <- matrix(nrow = 0, ncol = 2)
  colnames(random_genes) <- c("iteration", "average_correlation")
  for (i in (1:r)){ # Do this 100 times (could also try doing 1000 times if it doesn't take too long for the loop to run)
    genes <- sample(gene_names, k) # random sample from TST transcriptome
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

## D. Output file
write.csv(random_distributions_100, "../random_distrib_4-570size_100its_tpm_TST-2022_avg-correl.csv")  # change file name as appropriate (inc. iterations and range)
#random_distributions_100 <- read.csv("random_distrib_4-600size_100its_tpm_TST_avg-correl.csv")

## F. Plotting average cluster expression vs. frequency distribution
library(ggplot2)

ipa_avg_corr <- ipa_sep_reduced
ipa_avg_corr$avg_corr <- as.numeric(ipa_avg_corr$avg_corr)
random_distributions_100_df <- as.data.frame(random_distributions_100) #plotting only works with data frame

# load files if not running steps A-C above
#ipa_avg_corr <- read.csv("tpm_TST_avg-correl.csv") # output file from Script 1
#random_distributions_100 <- read.csv("random_distrib_4-600size_100its_tpm_TST_avg-correl.csv") # read in output file from Script 2 if not running steps A-D above
#random_distributions_100_df <- as.data.frame(random_distributions_100) #plotting only works with data frame

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

# Plotting the distributions only
random_distributions_plot <- ggplot(random_distributions_100_df, aes(X97.72, size)) + 
  geom_point(color = "palegreen3") +
  xlab("Average correlation coefficient") + ylab("Network size (number of target genes)") + #replace axes names here
  #ggtitle("Average Cluster Correlation Coefficient") + # replace title here
  geom_point(data = random_distributions_100_df, aes(X84.13, size), color = "steelblue2") +
  geom_point(data = random_distributions_100_df, aes(X99.86, size), color = "sienna2")+
theme_light(base_size = 14)
plot(random_distributions_plot)


