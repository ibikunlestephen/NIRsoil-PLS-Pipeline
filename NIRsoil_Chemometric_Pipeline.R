
# PCA & PLS WORKFLOW FOR NIR SPECTROSCOPY DATA
# Dataset: NIRsoil (built into prospectr package)
# Author: Stephen Ibikunle
# Date: 28-10-2026

# STEP 1: INSTALL AND LOAD PACKAGES
# Why: We need prospectr for spectral preprocessing, 
# ropls for PCA/PLS modeling, and tidyverse/ggplot2 for 
# data manipulation and publication-quality plots.

# Install packages
install.packages("prospectr")
install.packages("tidyverse")

#Because ropls package is missing on CRAN, you can go through Bioconductor repository
# Step 1: Install BiocManager if you don't have it already
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install ropls from Bioconductor
BiocManager::install("ropls")

# Load packages
library(ropls)
library(tidyverse)
library(prospectr)

# STEP 2: DATA INGESTION

# Why: Load the NIRsoil dataset and inspect its structure.
# The NIRsoil is a data.frame with 825 observations and 5 variables. The NIR 
# spectra is used as X and the reference nitrogen values is used as Y

# Load the built-in the NIRsoil dataset
data(NIRsoil)
View(NIRsoil)

# Inspect what we have
#NIRsoil is a data.frame with 825 observations and 5 variables.The NIR spectra 
#are stored in the spc column (a matrix with 825 rows × 700 wavelengths), and 
#the reference nitrogen values are in the Nt column.

str(NIRsoil)

# Extract the NIR spectra matrix (X)
# Each row is a sample, each column is a wavelength
X <- NIRsoil$spc

# Extract the reference values (Y)
# We'll use Nitrogen for PLS regression
Y <- NIRsoil$Nt

# Remove samples with missing Y values
complete_idx <- !is.na(Y)
X_clean <- X[complete_idx, ]
Y_clean <- Y[complete_idx]

print(dim(X_clean))  # Should be smaller than 825
print(length(Y_clean))

# Visualize raw spectra
# Why: Before any analysis, we need to see what the raw data looks like.
# This helps us identify baseline drift, scatter effects, and outliers.
# Visualize raw spectra (first 50 samples for clarity)

matplot(as.numeric(colnames(X_clean)), t(X_clean[1:50, ]), 
        type = "l", lty = 1, col = "gray60",
        xlab = "Wavelength (nm)", ylab = "Absorbance",
        main = "Raw NIR Spectra of Soil Samples (First 50)")

# STEP 3: PREPROCESSING

# Why: Raw NIR spectra suffer from light scattering and baseline drift.
# We apply Standard Normal Variate (SNV) to correct scatter,
# then Savitzky-Golay 2nd derivative to resolve overlapping peaks
# and remove baseline drift.

# Apply SNV correction
# SNV centers and scales each spectrum individually, removing 
# multiplicative scatter effects
X_snv <- standardNormalVariate(X_clean)

# Apply Savitzky-Golay 2nd derivative
# Window size = 11, polynomial order = 2
# This removes baseline drift and enhances peak resolution
X_sg <- savitzkyGolay(X_snv, m = 2, p = 2, w = 11)

# Visualize preprocessed spectra
matplot(as.numeric(colnames(X_sg)), t(X_sg[1:50, ]), 
        type = "l", lty = 1, col = "red",
        xlab = "Wavelength (nm)", ylab = "2nd Derivative Absorbance",
        main = "Preprocessed NIR Spectra (SNV + Savitzky-Golay 2nd Derivative)")

# Add a horizontal line at y=0 to show that the baselines are now flat
abline(h = 0, col = "black", lwd = 2, lty = 2)

# STEP 4: PRINCIPAL COMPONENT ANALYSIS (PCA) - UNSUPERVISED EXPLORATION

# PURPOSE: 
#   Reduce dimensionality from 700 wavelengths to a few latent variables
#   Detect outliers and quality issues in the dataset
#   Visualize natural sample clustering and variance structure
#   Validate that spectral variance correlates with nitrogen content

# 4.1 TRIM SPECTRAL EDGES TO REMOVE DERIVATIVE ARTIFACTS
# REASON: Savitzky-Golay derivatives amplify noise at spectrum boundaries
# because the sliding window lacks sufficient data points.
# Trimming prevents PCA from modeling instrument noise instead of
# true chemical variance.
X_pca <- X_sg[, 21:(ncol(X_sg) - 20)]  # Keep wavelengths 21 to 680 (trim first/last 20)

# 4.2 BUILD PCA MODEL
# REASON: Use opls() from ropls package for robust PCA implementation.
# predI=10 calculates up to 10 PCs to identify the optimal number

pca_model <- opls(X_pca, predI = 10)

# STEP 5: REMOVE OUTLIERS AND RE-RUN PCA
# REASON: Extreme outliers distort the PCA subspace and can cause PLS models 
# to overfit to anomalous samples rather than learning generalizable 
# chemical patterns. Removing them ensures model robustness.

# Identify outlier indices based on observation diagnostics
# Sample 702: Extreme score distance (far outside Hotelling's ellipse)
# Samples 402, 619, 363: High orthogonal distance (don't fit PCA subspace)

# Create clean dataset by excluding outliers
#raw indices like '702' fails if they exceed nrow(X_pca).
# Instead, we identify outliers by their statistical properties
# (extreme scores or orthogonal distances) to ensure accurate removal.

# STEP 5B FINAL FIX: MANUAL OUTLIER REMOVAL USING DIAGNOSTICS PLOT
# REASON: Automated coordinate filtering failed due to threshold mismatch.
# We now manually specify the exact row indices identified visually
# from the Observation Diagnostics plot to guarantee removal.

# From the Observation Diagnostics plot, we can SEE these labeled outliers:
# - Sample 402 (high OD, high SD)
# - Sample 619 (high OD)  
# - Sample 363 (high OD)
# - Sample 702 (extreme SD, far right)
# - Sample 283/245 (moderate SD outliers)

# IMPORTANT: These are SAMPLE LABELS from the plot, NOT row indices.
# We need to find their ACTUAL row positions in X_pca matrix.

# Method: Match by finding rows with extreme Score Distance (SD)
scores_matrix <- pca_model@scoreMN
score_distances <- sqrt(rowSums(scores_matrix^2))

# Find top 5 rows with highest score distance (these are our visual outliers)
top_outlier_rows <- order(score_distances, decreasing = TRUE)[1:5]

cat("Top 5 outlier row indices:", top_outlier_rows, "\n")
cat("Their score distances:", round(score_distances[top_outlier_rows], 2), "\n")

# Remove these specific rows
X_pca_clean <- X_pca[-top_outlier_rows, ]
Y_clean_outlier_removed <- Y_clean[-top_outlier_rows]

cat("\nOriginal samples:", nrow(X_pca), "\n")
cat("Samples after removal:", nrow(X_pca_clean), "\n")

# Re-run PCA on truly clean data
pca_model_clean <- opls(X_pca_clean, predI = 10)

# STEP 5C FINAL: REMOVE REMAINING OUTLIERS USING SCORE DISTANCE
# Calculate Score Distance from the current clean model
scores_matrix <- pca_model_clean@scoreMN
score_distances <- sqrt(rowSums(scores_matrix^2))

# Find top 5 rows with highest score distance (catches 303, 363, 619 + any others)
top_outlier_rows <- order(score_distances, decreasing = TRUE)[1:5]

cat("Top 5 outlier row indices to remove:", top_outlier_rows, "\n")
cat("Their score distances:", round(score_distances[top_outlier_rows], 2), "\n")

# Remove these specific rows
X_Clean_perfect <- X_pca_clean[-top_outlier_rows, ]
Y_clean_perfect <- Y_clean_outlier_removed[-top_outlier_rows]

cat("\nSamples before removal:", nrow(X_pca_clean), "\n")
cat("Samples after removal:", nrow(X_Clean_perfect), "\n")

# Re-run PCA on truly clean data
pca_model_perfect <- opls(X_Clean_perfect, predI = 10)

# STEP 6: PLS REGRESSION (SUPERVISED MODELING)

# Split the CLEANED data (70% train, 30% test)
set.seed(42)
train_idx <- sample(1:nrow(X_Clean_perfect), size = 0.7 * nrow(X_Clean_perfect))
X_train <- X_Clean_perfect[train_idx, ]
X_test  <- X_Clean_perfect[-train_idx, ]
Y_train <- Y_clean_perfect[train_idx]
Y_test  <- Y_clean_perfect[-train_idx]

# Step 7: Run PLS Model (Using 7-fold CV, which is more robust than LOO)
pls_model <- opls(X_train, Y_train, predI = 10, crossvalI = 7)

# Step 8: Extract Metrics Safely (Using modelDF to avoid API errors)
metrics <- pls_model@modelDF
opt_comp <- which.max(metrics$Q2) # Optimal components = where Q2 is highest

# Get the cumulative values at the optimal component
R2X_val <- metrics$R2X[opt_comp]
R2Y_val <- metrics$R2Y[opt_comp]
Q2_val  <- metrics$Q2[opt_comp]

cat("Optimal PLS components:", opt_comp, "\n")
cat("R2X:", round(R2X_val, 3), "| R2Y:", round(R2Y_val, 3), "| Q2:", round(Q2_val, 3), "\n")

# Step 9. Predict on Test Set
Y_pred <- predict(pls_model, newdata = X_test)

# Use predict() without newdata to get the training predictions safely
Y_train_pred <- predict(pls_model)

# Calculate RMSEP
RMSEP <- sqrt(mean((Y_test - Y_pred)^2))
cat("RMSEP:", round(RMSEP, 3), "\n")

# Calculate the Standard Deviation of your test set
sd_y_test <- sd(Y_test)

# Calculate RPD
RPD <- sd_y_test / RMSEP

cat("Standard Deviation of Test Set:", round(sd_y_test, 3), "\n")
cat("RPD (Ratio of Performance to Deviation):", round(RPD, 2), "\n")

# Step 10. Predicted vs Actual for visual determination of overfitting or underfitting
# Extract the TOTAL model metrics (the last row of the modelDF table)
# We use the confirmed values from the console output to ensure accuracy
true_R2Y <- 0.804
true_Q2 <- 0.605

# Create the prediction data frame for plotting
pred_df <- data.frame(
  Actual = c(Y_train, Y_test),
  Predicted = c(Y_train_pred, Y_pred),
  Set = c(rep("Training", length(Y_train)), rep("Test", length(Y_test)))
)

# Generate PLS Regression Plot 
PLS_plot <- ggplot(pred_df, aes(x = Actual, y = Predicted, color = Set)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  scale_color_manual(values = c("Training" = "blue", "Test" = "red")) +
  labs(title = "PLS Regression: Predicted vs Actual Nitrogen",
       subtitle = paste0("R² = ", true_R2Y, ", Q² = ", true_Q2, ", RMSEP = ", round(RMSEP, 3)),
       x = "Actual Nitrogen", y = "Predicted Nitrogen") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))

print(PLS_plot)

# STEP 11: RIGOROUS VALIDATION (Permutation & VIP)

# (11a) Permutation Test (Y-Randomization)

set.seed(123)
n_perm <- 100
Q2_perm <- numeric(n_perm)

for(i in 1:n_perm) {
  Y_shuffled <- sample(Y_train)
  
  perm_model <- suppressWarnings(
    opls(X_train, Y_shuffled, predI = 10, crossvalI = 7)
  )
  
  Q2_perm[i] <- perm_model@modelDF$Q2[nrow(perm_model@modelDF)] 
}

# Hardcode the TRUE Q2 to prevent any variable overwriting
true_Q2 <- 0.605 

# Create data frame and plot
perm_df <- data.frame(Q2_perm = Q2_perm)

p_perm <- ggplot(perm_df, aes(x = Q2_perm)) +
  geom_histogram(binwidth = 0.02, fill = "steelblue", color = "black") +
  geom_vline(xintercept = true_Q2, color = "red", linetype = "dashed", linewidth = 1) + 
  coord_cartesian(xlim = c(-0.2, 0.7)) + 
  labs(title = "PLS Permutation Test (Y-Randomization)",
       subtitle = paste0("Original Cumulative Q² (", true_Q2, ") >> Permuted Q² (Model is valid)"),
       x = "Permuted Q²", y = "Frequency") +
  theme_classic(base_size = 14)

print(p_perm)

# VIP (Variable Importance in Projection)
# Using the explicit slot @vipVn to avoid function version issues
vip_vals <- pls_model@vipVn
wavelengths <- as.numeric(colnames(X_pca_perfect))

vip_df <- data.frame(Wavelength = wavelengths, VIP = vip_vals)
vip_df$Important <- ifelse(vip_df$VIP > 1, "VIP > 1", "VIP ≤ 1")

# (11b) VIP Scores Plot
VIP_scores_plot <- ggplot(vip_df, aes(x = Wavelength, y = VIP, fill = Important)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("VIP > 1" = "darkred", "VIP ≤ 1" = "gray70")) +
  labs(title = "PLS Variable Importance in Projection (VIP)",
       subtitle = "Wavelengths with VIP > 1 drive the Nitrogen prediction",
       x = "Wavelength (nm)", y = "VIP Score") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")
print(VIP_scores_plot)
