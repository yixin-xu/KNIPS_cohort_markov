require(data.table)

generate_net_benefit_lapply <- function(input_parameters, lambda = 20000) {
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
  
  
  # Implant costs (transpose to keep convention of n_implants, n_samples)
  implant_costs <- t(input_parameters[, grepl("implant_cost", colnames(input_parameters))])
  rownames(implant_costs) <- treatment_names
  # Build an array to store the cohort vector at each cycle
  # Each cycle, implant and PSA sample has a single cohort vector of length n_states
  cohort_vectors <- array(0, dim = c(n_cycles, n_treatments, n_samples,  n_states), 
                          dimnames = list(NULL, treatment_names, NULL, state_names))
  
  mortality <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "mortality")
  primary_mortality=rep(abs(rnorm(n_samples, mean = as.numeric(mortality[which(grepl(paste0(initial_age," ", gender),mortality$Primary)),"estimate...3"]), 
                                  sd = (as.numeric(mortality[which(grepl(paste0(initial_age," ", gender),mortality$Primary)),"95%CI high...6"])-as.numeric(mortality[which(grepl(paste0(initial_age," ", gender),mortality$Primary)),"95%CI low...5"])/2*1.96))), each =n_treatments)
  
  # Assume everyone starts in the post_thr state
  cohort_vectors[1, , , "State Post TKR <3 years"] <- 1-primary_mortality
  # All other proportions start at zero
  cohort_vectors[1, , , "State Death"] <- primary_mortality
  
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
  
  # Pre-calculate the discount vector to reduce runtime
  discount_vector <- (1 / 1.035)^rep(c(0:(n_cycles-1)), each = 1)
  

  # Main model code
  # Use lapply instead of loop
  lapply(c(1:n_treatments), function(i_treatment){
    # Pre-index to reduce runtime
    transition_matrices_tr <- transition_matrices[, i_treatment, , , ]
    cohort_vectors_tr <- cohort_vectors[, i_treatment, , ]
    cycle_costs_tr <- cycle_costs[, i_treatment, ]
    cycle_qalys_tr <- cycle_qalys[, i_treatment, ]
    treatment_costs_tr <- implant_costs[i_treatment, ]
    total_costs_tr <- total_costs[i_treatment, ]
    total_qalys_tr <- total_qalys[i_treatment, ]
    state_costs_tr <- state_costs[ ,i_treatment, ]
    # In this case state qalys are the same for all treatments/implants 
    # but in general allow this for optimization
    state_qalys_tr <- state_qalys[,i_treatment, ]

    
    # Loop over the cycles
    # Cycle 1 is already defined so only need to update cycles 2:n_cycles
    
    # Loop over the PSA samples
    for (i_sample in 1:n_samples) {
      
      transition_matrices_tr_sample <- transition_matrices_tr[, i_sample, , ]
      
      cohort_vectors_tr_sample <- cohort_vectors_tr[, i_sample, ]
      
      # Loop over the cycles
      # Cycle 1 is already defined so only need to update cycles 2:n_cycles
      for(i_cycle in 2:n_cycles){
        # Markov update
        # Multiply previous cycle's cohort vector by transition matrix
        # i_e_ pi_j = pi_(j-1)*P
         cohort_vectors_tr_sample[i_cycle, ] <- 
          cohort_vectors_tr_sample[i_cycle - 1, ] %*% transition_matrices_tr_sample[i_cycle - 1, , ]
      } 
      # Use cohort vectors to calculate cycle costs and qalys
      cycle_costs_tr[, i_sample] <- cohort_vectors_tr_sample %*% state_costs_tr[i_sample, ]
      cycle_qalys_tr[, i_sample] <- cohort_vectors_tr_sample %*% state_qalys_tr[i_sample, ]
      # Sum  and discount to get total costs and qalys
      # Add implant costs to total costs
      total_costs_tr[i_sample] <- treatment_costs_tr[i_sample] + cycle_costs_tr[, i_sample] %*% discount_vector
      total_qalys_tr[i_sample] <- cycle_qalys_tr[, i_sample] %*% discount_vector
    }
    
    
    
    return(list(total_qalys = total_qalys_tr, total_costs = total_costs_tr))
    
  }) -> output_list # End lapply
  
  
  names(output_list) <- treatment_names 
  
  # Reconvert result to a matrix
  total_costs <- sapply(treatment_names, function(treatment_names) {total_costs[treatment_names, ] <- output_list[[treatment_names]]$total_costs})
  total_qalys <- sapply(treatment_names, function(treatment_names) {total_qalys[treatment_names, ] <- output_list[[treatment_names]]$total_qalys})
  # sapply inverts the matrices to uninvert
  total_costs <- t(total_costs)
  total_qalys <- t(total_qalys)
  
  # Calculate net benefit and incremental net benefit at
  # willingness-to-pay that was supplied
  net_benefit <- total_qalys * lambda - total_costs
  incremental_net_benefit <- net_benefit - net_benefit[1, ]
  
  return(list("total_costs" = total_costs,
              "total_qalys" = total_qalys,
              "net_benefit" = net_benefit,
              "incremental_net_benefit" = incremental_net_benefit))
}