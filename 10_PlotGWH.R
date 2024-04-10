#### Script to plot all Genome-Wide heterozygosity analysis
#### sliding windows GWH, angsd global GWH, statistical analysis
#### cetaceans GWH

setwd("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Heterozygosity")
library(ggplot2)
library(dplyr)
library(readxl)

#### Genome-wide heterozygosity in sliding windows #### 

#hector1mb <- read.table("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Heterozygosity/hectors_1mb.windowed.pi", header = TRUE)

#maui1mb <- read.table("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/MauiFull/Heterozygosity/maui_1mb.windowed.pi", header = TRUE)

maui10kb <- read.table("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/MauiFull/Heterozygosity/maui_10kb.windowed.pi", header = TRUE)

hectors10kb <- read.table("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/Heterozygosity/hectors_10kb.windowed.pi", header = TRUE)

maverage_pi <- mean(maui10kb$PI)

haverage_pi <- mean(hectors10kb$PI)

color_vector <- rep(c("darkblue", "lightblue"), length.out = length(unique(maui10kb$CHROM)))

chromosomes <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21")

plotm <- ggplot(maui10kb, aes(x = 1:nrow(maui10kb), y = PI, color = CHROM)) +
  geom_line(size = 0.8) +
  scale_x_continuous(breaks = NULL, labels = maui10kb$CHROM) +
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(y = "Heterozygosity (π)", x = NULL, title = "Māui") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_manual(values = c("coral1", "red", "coral1", "red", "coral1", "red", "coral1", "red", "coral1", "red"
                                , "coral1", "red", "red", "coral1", "coral1", "red", "coral1", "red", "coral1", "red", "coral1"))
ggsave("mauired.png", plotm, width = 8, height = 6, dpi = 300)

ploth <- ggplot(hectors10kb, aes(x = 1:nrow(hectors10kb), y = PI, color = CHROM)) +
  geom_line(size = 0.8) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(y = "Heterozygosity (π)", x = NULL, title = "Hector's") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_manual(values = c("cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3"
                                , "cornflowerblue", "blue3", "blue3", "cornflowerblue", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue"))
ggsave("hectors.png", ploth, width = 8, height = 6, dpi = 300)


plot_hist <- ggplot(maui10kb, aes(x = PI)) +
  geom_histogram(binwidth = 0.001, fill = "blue", color = "black") +
  labs(x = "Heterozygosity (π)", y = "Count", title = "Histogram of Per-Window Heterozygosity Levels") +
  theme_minimal()

# Save the histogram as an image
ggsave("mhistogram.png", plot_hist, width = 8, height = 6, dpi = 300)

hplot_hist <- ggplot(hectors10kb, aes(x = PI)) +
  geom_histogram(binwidth = 0.001, fill = "blue", color = "black") +
  labs(x = "Heterozygosity (π)", y = "Count") +
  theme_minimal()

ggsave("hechistogram.png", plot_hist, width = 8, height = 6, dpi = 300)

#### AUTOCORRELATION in Hector's and Māui genomes ####

# Autocorrelation plot for Maui's data
acf_maui <- acf(maui10kb$PI, lag.max = 100, plot = FALSE)
plot_acf_maui <- ggplot(data.frame(lag = acf_maui$lag, acf = acf_maui$acf), aes(lag, acf)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  labs(y = "Autocorrelation", x = "Lag", title = "Autocorrelation for Maui's Data") +
  theme_minimal()

# Autocorrelation plot for Hector's data
acf_hectors <- acf(hectors10kb$PI, lag.max = 100, plot = FALSE)
plot_acf_hectors <- ggplot(data.frame(lag = acf_hectors$lag, acf = acf_hectors$acf), aes(lag, acf)) +
  geom_bar(stat = "identity", width = 0.5, fill = "cornflowerblue") +
  labs(y = "Autocorrelation", x = "Lag", title = "Autocorrelation for Hector's Data") +
  theme_minimal()

#### Angsd heterozygosity ####
#Hector's angsd genome wide heterozygosity
a<-scan("hest.ml")
a[2]/sum(a)

#Māui angsd genome wide heterozygosity
b<-scan("mest.ml")
b[2]/sum(b)


#### T-test of Hector's Vs Māui mean genome wide heterozygosity #### 

# Combine data from both individuals into a single dataset
maui10kb$Individual <- "Māui"
hectors10kb$Individual <- "Hector's"
combined_data <- rbind(maui10kb, hectors10kb)

# Perform a t-test on the merged dataset
t_test_results <- t.test(combined_data$PI ~ combined_data$Individual)

# Print the t-test results
print(t_test_results)

# Extract and print the p-value
p_value <- t_test_results$p.value
cat("P-Value:", p_value, "\n")

# Check if the difference is statistically significant (e.g., p-value < 0.05)
if (p_value < 0.05) {
  cat("The difference in heterozygosity is statistically significant.\n")
} else {
  cat("There is no statistically significant difference in heterozygosity.\n")
}

#### Mann Whitney U-test Hector's Vs Māui mean genome wide heterozygosity #### 
# Load the required libraries if not already loaded
# install.packages("dplyr") # Install 'dplyr' if not already installed
library(dplyr)

# Create the combined dataset
maui10kb$Individual <- "Māui"
hectors10kb$Individual <- "Hector's"
combined_data <- rbind(maui10kb, hectors10kb)

# Perform a Mann-Whitney U test on the merged dataset
mann_whitney_results <- wilcox.test(PI ~ Individual, data = combined_data)

# Print the Mann-Whitney U test results
print(mann_whitney_results)

# Extract and print the p-value
p_value <- mann_whitney_results$p.value
cat("P-Value:", p_value, "\n")

# Check if the difference is statistically significant (e.g., p-value < 0.05)
if (p_value < 0.05) {
  cat("The difference in heterozygosity is statistically significant.\n")
} else {
  cat("There is no statistically significant difference in heterozygosity.\n")
}

#### BOXPLOT of mean Heterozygosity Hector's Vs Māui #### 
library(ggplot2)

# Calculate genome-wide mean and standard deviation for each individual
maverage_pi <- 0.0007 #Value obtained from angsd analysis
haverage_pi <- 0.0011 #Value obtained from angsd analysis
mstd_deviation_pi <- sd(maui10kb$PI)
hstd_deviation_pi <- sd(hectors10kb$PI)

# Create a data frame for mean and standard deviation
summary_stats <- data.frame(
  Individual = c("Māui", "Hector's"),
  Mean = c(maverage_pi, haverage_pi),
  StdDeviation = c(mstd_deviation_pi, hstd_deviation_pi)
)

# Create a vertical box plot with standard deviation whiskers
box_plot <- ggplot(combined_data, aes(x = Individual, y = PI)) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.0015)) +
  labs(x = NULL,y = "Heterozygosity (π)",title = "Box Plot of Mean Heterozygosity Comparison between Māui and Hector's") +
  theme_minimal()

library(ggplot2)

# Create separate data subsets for Hector's and Maui individuals
hectors_data <- subset(combined_data, Individual == "Hector's")
maui_data <- subset(combined_data, Individual == "Maui")

# Create the boxplot with outliers removed and custom colors
plot <- ggplot() +
  geom_boxplot(data = hectors_data, aes(x = Individual, y = PI), fill = "blue") +
  geom_boxplot(data = maui_data, aes(x = Individual, y = PI), fill = "red") +
  scale_y_continuous(limits = c(0, 0.0015)) +
  labs(x = NULL, y = "Heterozygosity (π)", title = "Box Plot of Mean Heterozygosity Comparison between Māui and Hector's") +
  theme_minimal()

# Remove outliers
plot + coord_cartesian(ylim = c(0, 0.0015))

# Save the box plot as an image
ggsave("mean_std_boxplot.png", box_plot, width = 8, height = 6, dpi = 300)

#### Violin Plot of Heterozygosity Distribution between Māui and Hector's ####

violin_plot <- ggplot(combined_data, aes(x = Individual, y = PI, fill = Individual)) +
  geom_violin(trim = FALSE, adjust = 2) +  # Adjust the bandwidth parameter (try different values)
  scale_fill_manual(values = c("blue", "red")) +
  labs(x = NULL, y = "Heterozygosity (π)", title = "Violin Plot of Heterozygosity Distribution between Māui and Hector's") +
  theme_minimal() +
  ylim(0, 0.003)
ggsave("violin_plot.png", violin_plot, width = 8, height = 6, dpi = 300)

#### Genome heterozygosity from 10x Chromium linked reads #### 

maui10x <- read.table("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/10X_PSMC/mauixmaui_10kb.windowed.pi", header = TRUE)

plotm <- ggplot(maui10x, aes(x = 1:nrow(maui10x), y = PI, color = CHROM)) +
  geom_line(size = 0.8) +
  scale_x_continuous(breaks = NULL, labels = maui10x$CHROM) +
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(y = "Heterozygosity (π)", x = NULL, title = "Māui") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
  #scale_color_manual(values = c("cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3"
                                #, "cornflowerblue", "blue3", "blue3", "cornflowerblue", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue"))
ggsave("maui10x.png", plotm, width = 8, height = 6, dpi = 300)

hectors10x <- read.table("/nesi/nobackup/uoo02423/Sebastian/MergedAssemblies/Hectors42x/10X_PSMC/2hectors10x_10kb.windowed.pi", header = TRUE)

plotm <- ggplot(hectors10x, aes(x = 1:nrow(hectors10x), y = PI, color = CHROM)) +
  geom_line(size = 0.8) +
  scale_x_continuous(breaks = NULL, labels = hectors10x$CHROM) +
  scale_y_continuous(limits = c(0, 0.1)) +
  labs(y = "Heterozygosity (π)", x = NULL, title = "hectors") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_color_manual(values = c("cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3"
                                , "cornflowerblue", "blue3", "blue3", "cornflowerblue", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue", "blue3", "cornflowerblue"))
ggsave("Hectors10x.png", plotm, width = 8, height = 6, dpi = 300)


#### Plot GWH of all available data from cetacean species #### 
# Load data from Excel
datacet <- read_excel("~/ResultsTables/GWHforGraph.xlsx")

# Order the data by heterozygosity (highest to lowest)

# Create the plot
pcet <- ggplot(datacet, aes(x = reorder(Species, Heterozygosity, decreasing = TRUE), y = Heterozygosity, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(title = "Genome-wide Heterozygosity by Species",
       x = "Species",
       y = "Genome-wide Heterozygosity") +
  scale_fill_manual(values = c("Hector's dolphin" = "blue",
                               "Māui dolphin" = "red",
                               "Other Species" = "gray")) +  # Set custom fill colors
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave("CetaceansGWH.png", plotm, width = 8, height = 6, dpi = 300)

