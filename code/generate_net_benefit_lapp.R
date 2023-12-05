require(reshape)
require(dplyr)

generate_net_benefit_df <- function(input_parameters, lambda = 20000) {
  # First generate components needed for simulation
  transition_matrices <- generate_transition_matrices(input_parameters, 
                                                      treatment_names = treatment_names, 
                                                      state_names = state_names,
                                                      initial_age = initial_age,
                                                      ini_age = ini_age,
                                                      final_age = final_age,
                                                      starting_age = starting_age,
                                                      gender = gender,
                                                      sensitivity = NULL)
  state_costs <- generate_state_costs(input_parameters,
                                      treatment_names = treatment_names, 
                                      state_names = state_names,
                                      initial_age = initial_age,
                                      final_age = final_age,
                                      starting_age = starting_age,
                                      gender = gender)
  state_qalys <- generate_state_qalys(input_parameters,
                                      treatment_names = treatment_names, 
                                      state_names = state_names,
                                      initial_age = initial_age,
                                      final_age = final_age,
                                      starting_age = starting_age,
                                      gender = gender)
  
  # Convert transition matrices array to a dataframe
  # Must be in cycle, implant, sample, "from state", "to states" order
  # Also remove dimnames so data.frame can index using numbers (important if going to C/C++)
  transition_matrices_temp <- transition_matrices
  dimnames(transition_matrices_temp) <- list(NULL, NULL, NULL, NULL, state_names)
  # List of data frames for transition from each state
  transition_matrices_temp_df <- list()
  for(i_state in 1:n_states) {
    # Each stores the transition probabilities from i_state
    transition_matrices_temp_df[[i_state]] <- transition_matrices_temp[, , , , i_state]
    # Convert the multidimensional array to a data frame
    transition_matrices_temp_df[[i_state]] <- 
      melt(transition_matrices_temp_df[[i_state]], varnames = c("cycle", "implant", "sample", "from"))
    # Name the state to which you're transiting
    colnames(transition_matrices_temp_df[[i_state]])[5] <- state_names[i_state]
  }
  # Combine the data frames
  transition_matrices_df <- do.call(cbind, transition_matrices_temp_df)
  # Only keep unique columns
  transition_matrices_df <- transition_matrices_df[, unique(colnames(transition_matrices_df))]
  # Sort so that from asce
  transition_matrices_df %>% arrange(cycle, implant, sample, from)
  
  # Check if we've set them up correctly
  #transition_matrices_df[with(transition_matrices_df, cycle == 5 & implant == 3 & sample == 101), c(4:8)]
  #transition_matrices[5, 3, 101, , ]
  
  
  
  # Store the cohort vectors as a data frame with one row for each cycle, implant and sample
  cohort_vectors <- data.frame(cycle = rep(c(1:n_cycles), n_treatments * n_samples),
                               implant = rep(c(1:n_treatments), each = n_cycles, n_samples),
                               sample = rep(c(1:n_samples), each = n_cycles * n_treatments))
  # One column for each states
  for(i_state in 1:n_states) {cohort_vectors <- cbind(cohort_vectors, rep(0, n_cycles * n_treatments * n_samples))}
  colnames(cohort_vectors)[4:(3 + n_states)] <- state_names
  
  mortality <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "mortality")
  primary_mortality=rep(abs(rnorm(n_samples, mean = as.numeric(mortality[which(grepl(paste0(initial_age," ", gender),mortality$Primary)),"estimate...3"]), 
                                  sd = (as.numeric(mortality[which(grepl(paste0(initial_age," ", gender),mortality$Primary)),"95%CI high...6"])-as.numeric(mortality[which(grepl(paste0(initial_age," ", gender),mortality$Primary)),"95%CI low...5"])/2*1.96))), each =n_treatments)
  # Assume everyone starts in the post_thr state
  cohort_vectors[cohort_vectors$cycle == 1, "State Post TKR <3 years"] <- 1-primary_mortality
  cohort_vectors[cohort_vectors$cycle == 1, "State Death"] <- primary_mortality
  # All other proportions start at zero by default when setting up data frame
  
  
  # Implant costs (transpose to keep convention of n_implants, n_samples)
  implant_costs <- t(input_parameters[, grepl("implant_cost", colnames(input_parameters))])
  rownames(implant_costs) <- treatment_names
  
  # Build an array to store the costs and QALYs accrued per cycle
  # One for each cycle, implant, and PSA sample
  # These will be filled in below in the main model code
  # Then discounted and summed to contribute to total costs and total QALYs
  cycle_costs <- array(dim = c(n_cycles, n_treatments, n_samples), 
                       dimnames = list(NULL, treatment_names, NULL))
  cycle_qalys <- array(dim = c(n_cycles, n_treatments, n_samples), 
                       dimnames = list(NULL, treatment_names, NULL))
  
  # Build arrays to store the total costs and total QALYs
  # There is one for each treatment and each PSA sample
  # These are filled in below using cycle_costs, 
  # implant_costs, and cycle_qalys
  total_costs <- array(dim = c(n_treatments, n_samples), 
                       dimnames = list(treatment_names, NULL))
  total_qalys <- array(dim = c(n_treatments, n_samples), 
                       dimnames = list(treatment_names, NULL))
  
  
  # The remainder of the cohort_vectors will be filled in by Markov updating below
  discount_vector <- (1 / 1.035)^rep(c(0:(n_cycles-1)), each = 1)
  
  # Main model code
  # Loop over the implants
  
  for (i_implant in 1:n_treatments) {
    # Loop over the PSA samples
    for (i_sample in 1:n_samples) {
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for (i_cycle in 2:n_cycles) {
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i_e_ pi_j = pi_(j-1)*P
        
        # Extremely slow due to type conversions
        cohort_vectors[with(cohort_vectors, cycle == i_cycle & implant == i_implant & sample == i_sample), c(4:(3 + n_states))] <-
          as.data.frame(as.numeric(as.vector(cohort_vectors[with(cohort_vectors, cycle == i_cycle - 1 & implant == i_implant & sample == i_sample), c(4:(3 + n_states))])) %*%
                          as.matrix(transition_matrices_df[with(transition_matrices_df, cycle == i_cycle - 1 & implant == i_implant & sample == i_sample), c(5:(4 + n_states))]))
        
      }
    }
  }
  
  
  lapply(c(1:n_treatments), function(i_treatment){
    # Pre-index to reduce runtime
    #cohort_vectors_tr <- cohort_vectors[i_treatment ] 
    cycle_costs_tr <- cycle_costs[, i_treatment, ]
    cycle_qalys_tr <- cycle_qalys[, i_treatment, ]
    treatment_costs_tr <- implant_costs[i_treatment, ]
    total_costs_tr <- total_costs[i_treatment, ]
    total_qalys_tr <- total_qalys[i_treatment, ]
    state_costs_tr <- state_costs[ ,i_treatment, ]
    # In this case state qalys are the same for all treatments/implants 
    # but in general allow this for optimization
    state_qalys_tr <- state_qalys[,i_treatment, ]
    
    
    for(i_sample in 1:n_samples) {
      cohort_vectors_tr_sample <- cohort_vectors[c(0:(n_cycles-1)) * (n_treatments * n_samples) +
                                                   (i_treatment - 1) * n_samples + i_sample, c(4:(3+n_states))]
      # Use cohort vectors to calculate cycle costs and qalys
      cycle_costs_tr[, i_sample] <- as.matrix(cohort_vectors_tr_sample[,  ]) %*% state_costs_tr[i_sample, ]
      cycle_qalys_tr[, i_sample] <- as.matrix(cohort_vectors_tr_sample[,  ]) %*% state_qalys_tr[i_sample, ]
      # Sum  and discount to get total costs and qalys
      # Add implant costs to total costs
      total_costs_tr[i_sample] <- treatment_costs_tr[i_sample] + cycle_costs_tr[, i_sample] %*% discount_vector
      total_qalys_tr[i_sample] <- cycle_qalys_tr[, i_sample] %*% discount_vector #%*% multi_factor
    }
    
    
    
    return(list(total_qalys = total_qalys_tr, total_costs = total_costs_tr))
    
  }) -> output_list # End lapply
  
  
  names(output_list) <- treatment_names 
  
  # Reconvert result to a matrix
  total_costs <- sapply(treatment_names, function(treatment_name) {total_costs[treatment_name, ] <- output_list[[treatment_name]]$total_costs})
  total_qalys <- sapply(treatment_names, function(treatment_name) {total_qalys[treatment_name, ] <- output_list[[treatment_name]]$total_qalys})
  # sapply inverts the matrices to uninvert
  total_costs <- t(total_costs)
  total_qalys <- t(total_qalys)
  lambda = 20000
  # Calculate net benefit and incremental net benefit at
  # willingness-to-pay that was supplied
  net_benefit <- total_qalys * lambda - total_costs
  incremental_net_benefit <- net_benefit - matrix(rep(net_benefit[1, ], each = n_treatments), nrow = n_treatments)
  incremental_costs <- total_costs[,]-matrix(rep(total_costs[1, ], each = n_treatments), nrow = n_treatments)
  incremental_qalys <- total_qalys[,]-matrix(rep(total_qalys[1, ], each = n_treatments), nrow = n_treatments)
  ICER = rowMeans(incremental_costs)/rowMeans(incremental_qalys)
  
  return(list("total_costs" = total_costs,
              "total_qalys" = total_qalys,
              "net_benefit" = net_benefit,
              "incremental_net_benefit" = incremental_net_benefit,
              "ICER" = ICER))
}

