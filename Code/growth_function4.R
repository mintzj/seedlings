# Fixed function
# Update 4/25: allow a, b, c to be a vector, one for each tree.
simulate_tree_growth <- function(t, n, a0, a1, a, b, c, ht_initial, sd, msurv, 
                                 break_mean, break_prop, form = "logistic") {
  # Initialize data frame to store results
  data_all <- data.frame(Tree_ID = integer(), Time_Step = integer(), Height = numeric(), 
                         Tree_State = integer())
  
  surv <- matrix(NA, nrow = n, ncol = t)
  # Initialize tree states and heights
  tree_state <- rep(1, n)  # 1 represents alive, 2 represents dead
  ht2 <- rep(NA, n)  # Initialize height vector
  brk <- rep(0, n)
  gr <- rnorm(n, 0, sd = sd)  # Generate random growth rates for each tree
  
  for (i in 1:t) {  # Time point t
    pi <- numeric(n)  # Initialize the survival probability vector
    
    for (j in 1:n) {  # Tree j
      if (tree_state[j] == 1) {  # If the tree is alive
        if (i == 1) {  # All trees survive in the first time step
          surv[j, i] <- 1
          ht2[j] <- ht_initial[j]  # Set initial height for newly born trees
        } else {
          tmp <- exp(a0 + a1 * ht2[j])
          if(tmp > 1E10){
            tmp <- 1E10
          }
          pi[j] <- min(msurv, (tmp / (1 + tmp)))  # Use the smaller of the two survival rates
          surv[j, i] <- rbinom(1, size = 1, prob = pi[j])
          
          if (surv[j, i] == 0) {  # If the tree dies
            tree_state[j] <- 2  # Update the status of the tree to dead
            ht2[j] <- NA  # Mark the height of dead trees as NA
          } else {
            # If tree survives, 
            
            break_test <- rbinom(1, size = 1, prob = break_prop)
            
            
            
            # Original Growth:
            # mu_growth <- a / (1 + exp((b - ht2[j]) / c))  # Default, logistic.
            mu_growth <- a[j] / (1 + exp((b[j] - ht2[j]) / c[j]))  # Default, logistic.
            
            # Ricker form is going to take more work.
            if(form == "ricker"){
              mu_growth <- a[j] * ht2[j] * exp(-b[j] * ht2[j])
            }
            
            if(form == "gompertz"){
              mu_growth <- a[j] * exp(-b[j] * c[i]^ht2[j])
              # Can be rewritten in this form:
              # mu_growth <- exp(-a * exp(b * (x - c)))
            }
            
            
            ht2[j] <- ht2[j] + mu_growth + 
              gr[j] + break_test * break_mean  
            # If tree breaks, break by the average amount for species.
            
            brk[j] <- break_test
            
            if( (ht2[j] < 0) |  is.na(ht2[j]) ){
              ht2[j] = 0
            }
          }
        }
      } else {  # If the tree is dead
        surv[j, i] <- NA  # Set survival to NA for dead trees
        ht2[j] <- NA  # Ensure the height of dead trees remains as NA
      }
    } ## closes the j loop (tree)
    
    # Create the data frame for each time step
    tree_ID <- 1:n  
    height <- ht2
    tree_state_temp <- tree_state  # Temporary variable for tree state
    
    data <- data.frame(Tree_ID = tree_ID, Time_Step = rep(i, n), Height = height, 
                       Tree_State = tree_state_temp, Break = brk)
    data_all <- rbind(data_all, data)  # Append new data to the existing data frame
  } ## closes i loop (time point)
  return(data_all)
}

summarise_results <- function(result, maxh = 150){
  maxheights <- result %>% group_by(Tree_ID) %>% summarise(maxh = max(Height, na.rm = TRUE))
  
  num_to_ht_i <- sum(maxheights$maxh > maxh)
  num_needed_i <- n_seedling / num_to_ht_i
  
  if (any(result$Height >= maxh, na.rm = TRUE)) {
    greater_than <- result[which(result$Height >= maxh), ]
    min_time_to_ht_i <- min(greater_than$Time_Step)

    med_time_to_ht_i <- median(greater_than$Time_Step, na.rm = TRUE)

  } else {
    min_time_to_ht_i <- med_time_to_ht_i <- NA
  }
  
  return(c(num_to_ht_i, num_needed_i, min_time_to_ht_i, med_time_to_ht_i))
}
