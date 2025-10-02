# Load necessary libraries
library(compositions)  # For compositional data analysis
library(MASS)          # For multivariate normal distribution

# Define new parameters with their respective weights
parameters <- c(
  `Aluminum [µg/l Al]` = 0.13732495,
  `Ammonium [mg/l NH4]` = 0.08499825,
  `Arsenic [µg/l As]` = 0.00785124,
  `Cadmium [µg/l Cd]` = 0.03023614,
  `Calcium [mg/l Ca]` = 0.0037585,
  `Chlorides [mg/l Cl]` = 0.12006373,
  `Conductivity at 20°C [µS/cm]` = 0.10301697,
  `Copper [mg/l Cu]` = 0.05903513,
  `Fluorides [mg/l F]` = 0.05134627,
  `Hardness [°F]` = 0.00831767,
  `Iron [µg/l Fe]` = 0.02670118,
  `Lead [µg/l Pb]` = 0.02825393,
  `Manganese [µg/l Mn]` = 0.00192557,
  `Nitrates [mg/l NO3]` = 0.02986598,
  `pH` = 0.03950439,
  `Sodium [mg/l Na]` = 0.08998053,
  `Sulfates [mg/l SO4]` = 0.02081114,
  `Turbidity [NTU]` = 0.08237665,
  `Vanadium [µg/l V]` = 0.03977873,
  `Zinc [µg/l Zn]` = 0.03485305
)

# Normalize parameters to sum to 1
parameters <- parameters / sum(parameters)

# Convert parameters to a data frame with one row
mean_comp <- as.data.frame(t(parameters))

# Convert to acomp object
mean_comp <- acomp(mean_comp)

# Number of simulations (samples)
n_simulations <- 10000000  # Adjust as needed

# Desired coefficient of variation (10%)
target_cv <- 0.10

# Number of components
D <- ncol(mean_comp)
n_ilr <- D - 1

# Transform mean composition to ilr coordinates
mean_ilr <- ilr(mean_comp)

# Set variance in ilr space
ilr_variance_value <- (target_cv)^2  # Approximation for small variances

# Covariance matrix in ilr space
cov_ilr <- diag(ilr_variance_value, n_ilr)

# Generate samples in ilr space
set.seed(123)  # For reproducibility
ilr_samples <- mvrnorm(n = n_simulations, mu = as.numeric(mean_ilr), Sigma = cov_ilr)

# Transform back to the simplex
comp_samples <- ilrInv(ilr_samples)

# Convert to data frame
weights_df <- as.data.frame(comp_samples)
colnames(weights_df) <- colnames(mean_comp)

# Check that the samples sum to 1
sums <- rowSums(weights_df)
if (!all(abs(sums - 1) < 1e-10)) {
  stop("Some samples do not sum to 1.")
}

# Compute actual CVs
sd_values <- apply(weights_df, 2, sd)
mean_values <- apply(weights_df, 2, mean)
cv_values <- sd_values / mean_values

# Print CVs
print("Actual Coefficient of Variation for each component:")
print(cv_values)

# Save the weights to a CSV file
write.csv(weights_df, "FG_LAKE_WEIGHTS_10M.csv", row.names = FALSE)

# Plot the distributions of the sampled weights
pdf("weight_distributions_FG_LAKE.pdf")
par(mfrow = c(5, 4))  # Updated to accommodate 20 parameters
for (i in 1:ncol(weights_df)) {
  hist(
    weights_df[, i],
    main = colnames(weights_df)[i],
    xlab = "Weight",
    col = "skyblue",
    border = "white"
  )
}
dev.off()

print("Weight sampling complete. Check 'FG_WELLS_WEIGHTS_10M.csv' and 'weight_distributions_FG_WELLS.pdf' for results.")

