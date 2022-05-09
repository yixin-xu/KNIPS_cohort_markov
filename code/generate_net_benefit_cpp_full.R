# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert model inputs to ouputs using Markov model

library(dplyr)
# Uses the Rcpp package
require(Rcpp)
# Compiles the C++ file for the Markov loop
code_directory <- "C:/Users/yx18392/Desktop/knips/cohort markov model"
Rcpp::sourceCpp(file = paste0(code_directory, "/rcpp_loop_full.cpp"))


generate_net_benefit <- function(input_parameters, lambda = 20000) {
  # First generate components needed for simulation
  transition_matrices <- generate_transition_matrices(input_parameters)
  transition_matrices_df <- convert_transition_matrices_to_df(transition_matrices)
  state_costs <- generate_state_costs(input_parameters)
  state_qalys <- generate_state_qalys(input_parameters)
  
  state_costs[is.na(state_costs)] <- 0
  state_qalys[is.na(state_qalys)] <- 0
  
  data_directory <- "C:/Users/yx18392/Desktop/knips/data/rate"
  mortality <- read_excel(paste0(data_directory, "/cohort model inputs.xlsx"), sheet = "mortality")
  # Implant costs (transpose to keep convention of n_implants, n_samples)
  implant_costs <- t(input_parameters[, grepl("implant_cost", colnames(input_parameters))])
  rownames(implant_costs) <- treatment_names
  
  
  # Build an array to store the cohort vector at each cycle
  # Store the cohort vectors as a data frame with one row for each cycle, implant and sample
  cohort_vectors <- data.frame(cycle = rep(c(1:n_cycles), n_treatments * n_samples),
                               treatments = rep(c(1:n_treatments), each = n_cycles, n_samples),
                               sample = rep(c(1:n_samples), each = n_cycles * n_treatments))
  
  # One column for each states
  for(i_state in 1:n_states) {cohort_vectors <- cbind(cohort_vectors, rep(0, n_cycles * n_treatments * n_samples))}
  colnames(cohort_vectors)[4:(3 + n_states)] <- state_names
  
  # Sort 
  cohort_vectors <- cohort_vectors %>% arrange(cycle, treatments, sample)
  
  
  # Assume everyone starts in the post_thr state
  primary_mortality=rep(rnorm(n_samples, mean = as.numeric(mortality[2,3]), 
                              sd = (as.numeric(mortality[2,6])-as.numeric(mortality[2,5]))/2*1.96), each =n_treatments)
  cohort_vectors[cohort_vectors$cycle == 1, "State Post TKR <3 years"] <- 1-primary_mortality
  cohort_vectors[cohort_vectors$cycle == 1, "State Death"] <- primary_mortality
  # All other proportions start at zero by default when setting up data frame
  
  
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
  
  # Function currently requires conversion to matrices and returns a matrix, not a dataframe
  cohort_vectors <- 
    (rcpp_loop_full(cohort_vectors_in = as.matrix(cohort_vectors), 
                    transition_matrices = as.matrix(transition_matrices_df),
                    n_cycles = n_cycles, n_implants = n_treatments, n_samples = n_samples, n_states = n_states))
   
  
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
      cycle_costs_tr[, i_sample] <- cohort_vectors_tr_sample[,  ] %*% state_costs_tr[i_sample, ]
      cycle_qalys_tr[, i_sample] <- cohort_vectors_tr_sample[,  ] %*% state_qalys_tr[i_sample, ]
      # Sum  and discount to get total costs and qalys
      # Add implant costs to total costs
      total_costs_tr[i_sample] <- treatment_costs_tr[i_sample] + cycle_costs_tr[, i_sample] %*% discount_vector
      total_qalys_tr[i_sample] <- cycle_qalys_tr[, i_sample] %*% discount_vector
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

