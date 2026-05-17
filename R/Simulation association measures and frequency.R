# Set reproducibility so results are identical each run
set.seed(123)

# clear workspace
rm(list = ls())

# Large corpus size (controls expected value scaling)
N <- 1e6

# Number of simulated observations (kept small for visualisation clarity)
max_points <- 10000

# Pre-allocate vectors for efficiency (stores simulated values)
O_values <- numeric(max_points) # observed frequency-like variable
t_values <- numeric(max_points) # test statistic values
z_values <- numeric(max_points)
mu_values <- numeric(max_points)
mi_values <- numeric(max_points)
mi2_values <- numeric(max_points)
mi3_values <- numeric(max_points)
ll_values <- numeric(max_points)
dice_values <- numeric(max_points)
logdice_values <- numeric(max_points)
logratio_values <- numeric(max_points)
ms_values <- numeric(max_points)

# Main simulation loop
for (i in 1:max_points) {
  # Step 1: generate O (frequency)
  O_11 <- sample(1:5000, 1)

  # Step 2: generate R1 and C1
  # Constrained to be >= O_11
  R_1 <- sample(O_11:5000, 1)
  C_1 <- sample(O_11:5000, 1)

  # Step 3: calculate the remaining values of the observed table
  O_21 <- C_1 - O_11
  O_12 <- R_1 - O_11
  R_2 <- N - R_1
  C_2 <- N - C_1
  O_22 <- R_2 - O_21

  # Step 4: expected value under null model
  # Scaled down heavily by large corpus size N
  E_11 <- (R_1 * C_1) / N
  E_21 <- (R_2 * C_1) / N
  E_12 <- (R_1 * C_2) / N
  E_22 <- (R_2 * C_2) / N

  # Step 4: compute test statistic
  t <- (O_11 - E_11) / sqrt(O_11)
  # z-score
  z <- (O_11 - E_11) / sqrt(E_11)
  # MU
  mu <- O_11 / E_11
  # MI
  mi <- log2(O_11 / E_11)
  # MI2
  mi2 <- log2(O_11^2 / E_11)
  # MI2
  mi3 <- log2(O_11^3 / E_11)
  # LL
  ll <- 2 *
    (O_11 *
      log(O_11 / E_11) +
      O_21 * log(O_21 / E_21) +
      O_12 * log(O_12 / E_12) +
      O_22 * log(O_22 / E_22))
  # dice
  dice <- (2 * O_11) / (R_1 + C_1)
  # log dice
  logdice <- 14 + log2((2 * O_11) / (R_1 + C_1))
  # log ratio
  logratio <- log2((O_11 * R_2) / (O_21 * R_1))
  # minimum sensitivity
  ms <- min((O_11 / C_1) / O_11 / R_1)

  # Store results for plotting
  O_values[i] <- O_11
  t_values[i] <- t
  z_values[i] <- z
  mu_values[i] <- mu
  mi_values[i] <- mi
  mi2_values[i] <- mi2
  mi3_values[i] <- mi3
  ll_values[i] <- ll
  dice_values[i] <- dice
  logdice_values[i] <- logdice
  logratio_values[i] <- logratio
  ms_values[i] <- ms
}

measures <- list(
  T_score = t_values,
  Z_score = z_values,
  MU = mu_values,
  MI = mi_values,
  MI_2 = mi2_values,
  MI_3 = mi3_values,
  LogLikelihood = ll_values,
  Dice = dice_values,
  LogDice = logdice_values,
  LogRatio = logratio_values,
  MS = ms_values
)

# depending on the simulated data, some samples output `Inf` or `NaN` values
# lapply(lapply(measures, is.finite), sum)

# we set those non-numeric values to zero
measures_clean <- lapply(measures, function(v) {
  v[!is.finite(v)] <- 0
  v
})

# lapply(lapply(measures_clean, is.finite), sum)

# add observed frequency data to the list object
measures_clean$Frequency <- O_values

# normalize values between 0 and 1
normalize <- function(x) (x - min(x)) / (max(x) - min(x))

measures_clean_norm <- lapply(measures_clean, normalize)

# boxplot All Measures and Frequency
par(mfrow = c(1, 1), cex.main=1.2)
boxplot(
  measures_clean_norm,
  col = "steelblue",
  main = "Association Measures (normalized) and Observed Frequency",
  sub = "10,000 Random Simulations"
)

# Create scatterplots
# note - this is the version before cleanup
# Minimum Sensitivity is not included in this chart

measures <- measures[names(measures) != "MS"]


par(mfrow = c(2, 5), oma = c(0, 0, 3, 0))

lapply(names(measures), function(name) {
  plot(
    O_values,
    measures[[name]],
    main = name,
    xlab = "Observed Frequency",
    ylab = paste(name, "value"), # "value",
    pch = 16,
    col = "steelblue"
  )
})

mtext(
  "Association Measures by Observed Frequency",
  outer = TRUE,
  cex = 1.5
)

