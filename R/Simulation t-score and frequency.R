# Set reproducibility so results are identical each run
set.seed(123)
 
# Large corpus size (controls expected value scaling)
N <- 1e6
 
# Number of simulated observations (kept small for visualization clarity)
max_points <- 100
 
# Pre-allocate vectors for efficiency (stores simulated values)
O_values <- numeric(max_points)  # observed frequency-like variable
t_values <- numeric(max_points)  # test statistic values
 
# Set plotting layout (single plot window)
par(mfrow = c(1,1))
 
# Main simulation loop
for (i in 1:max_points) {
 
 # Step 1: generate O (frequency)
O <- sample(1:5000, 1)
 
 # Step 2: generate R1 and C1
 # Constrained to be >= O
R1 <- sample(O:5000, 1)
C1<- sample(O:5000, 1)
 
 # Step 3: expected value under null model
 # Scaled down heavily by large corpus size N
 E <- (R1 * C1) / N
 
 # Step 4: compute test statistic
 t <- (O - E) / sqrt(O)
 
 # Store results for plotting
 O_values[i] <- O
 t_values[i] <- t
 
 # Step 5: freeze-frame visualization
 # Only update plot every 10 iterations OR at final step
 if (i %% 10 == 0 || i == max_points) {
 
   plot(
    O_values[1:i], t_values[1:i],
    
     # axis limits fixed so growth is visually comparable across frames
     xlim = c(1, 5000),
     #ylim = c(0, 100),
    
     # point style
     pch = 16,
     col = rgb(0, 0.4, 0.8, 0.5),
    
     # dynamic title showing progression
     main = paste("Growth of O vs t | n =", i),
    
     # axis labels
     xlab = "O (observed frequency)",
     ylab = "t statistic"
   )
 }
}
 
