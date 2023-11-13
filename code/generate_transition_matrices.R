# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters to age-dependent transition matrices
# Output is a 3-dimensional array of n_samples x n_treatments x n_states x n_states


generate_transition_matrices <- function(input_parameters, 
                                         treatment_names = treatment_names, 
                                         state_names = state_names,
                                         initial_age = initial_age,
                                         ini_age = ini_age,
                                         final_age = final_age,
                                         starting_age = starting_age,
                                         gender = gender,
                                         n_cycles=50, sensitivity = NULL) {
  
  n_treatments <- length(treatment_names)
  n_samples <- dim(input_parameters)[1]
  n_states <- length(state_names)
  
  # create transition matrix 
  transition_matrices <- array(0, dim = c(n_cycles, n_treatments, n_samples, n_states, n_states),
                               dimnames = list(NULL, treatment_names, NULL, state_names, state_names))
  
  
  lifetime <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "UK_lifetables")
  mortality_90d_raw <- as.matrix(read_excel(paste0(data_directory,"/KNIPS Main input data.xlsx"), sheet = "90d_mortality"))
  mortality_row_index <- which(mortality_90d_raw[, "age_upper"] == final_age & mortality_90d_raw[, "gender"] == gender)
  
  # Mortality probability is sampled from a beta distribution
  # Each sample is converted to a rate
  tkr_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "tkr_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "tkr_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "tkr_prop_dead"]))
  first_mortality <- rbeta(n_samples, shape1 = tkr_mortality_90d_beta_params[1], shape2 = tkr_mortality_90d_beta_params[2])
  
  revision_mortality_90d_beta_params <- as.numeric(mortality_90d_raw[mortality_row_index, "revision_number"]) * c(
    as.numeric(mortality_90d_raw[mortality_row_index, "revision_prop_dead"]),
    1 - as.numeric(mortality_90d_raw[mortality_row_index, "revision_prop_dead"]))
  second_mortality <- rbeta(n_samples, shape1 = revision_mortality_90d_beta_params[1], shape2 = revision_mortality_90d_beta_params[2])
  higher_mortality <- second_mortality
   # HT: All the death probabilities should be the same. The i_cycle in the Markov loop already counts the increased age by the time patients enter the later states
  # For example, using age 62 for "State Post TKR >=3 years < 10 years" means that patients entering "State Post TKR <3 years" at age 60
  # will be 62 years old (i.e. i_cycle = 3) when they reach "State Post TKR >=3 years < 10 years" but their mortality will come from lifetime[65, ] since it's 62+i_cycle
  # We discussed this before so please try to fix. We can talk through it again in January if still not clear.
  # death rate
  
  
  for(i_cycle in 1:n_cycles) {
    for(treatment_name in treatment_names) {
      transition_matrices[i_cycle, treatment_name, , "State Post TKR <3 years", "State Early revision"] <- 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))*(1-first_mortality)
      transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State middle revision"] <- 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))*(1-first_mortality)
      transition_matrices[i_cycle, treatment_name, , "State Post TKR >=10 years", "State late revision"] <- 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))*(1-first_mortality)
      
      
      transition_matrices[i_cycle, treatment_name, , "State Early revision", "State second revision"] <- 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))*(1-second_mortality)
      transition_matrices[i_cycle, treatment_name, , "State middle revision", "State second revision"] <- 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))*(1-second_mortality)
      transition_matrices[i_cycle, treatment_name, , "State late revision", "State second revision"] <- 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))*(1-second_mortality)
      
      #if(!is.null(sensitivity)) {
        #if(sensitivity == "remove_2nd_higher_rate"){
          #transition_matrices[i_cycle, treatment_name, , "State Early revision", "State second revision"] = 
            #transition_matrices[i_cycle, treatment_name, , "State middle revision", "State second revision"]=
            #transition_matrices[i_cycle, treatment_name, , "State late revision", "State second revision"]= rep(NaN, each = n_samples)
       # }
     # }
          
      
      #transition_matrices[i_cycle, treatment_name, , "State second revision", "State second revision"] <- 
        #(1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))*(1-higher_mortality)
      
      # Annual probability of death the same for all states
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))*first_mortality
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))*first_mortality
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=10 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))*first_mortality
      
      transition_matrices[ i_cycle, treatment_name, , "State Early revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))*second_mortality
      transition_matrices[ i_cycle, treatment_name, , "State middle revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))*second_mortality
      transition_matrices[ i_cycle, treatment_name, , "State late revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))*second_mortality
      
      
      
      transition_matrices[ i_cycle, treatment_name, , "State second revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[starting_age + i_cycle - 1, gender])), n_samples) + 
        (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))*higher_mortality
      
      if(is.infinite(as.numeric(final_age))) {
        transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=10 years", "State Death"] = rep(1-exp(-as.numeric(lifetime[101, gender])), n_samples) + 
          (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))*first_mortality
        transition_matrices[ i_cycle, treatment_name, , "State late revision", "State Death"] = rep(1-exp(-as.numeric(lifetime[101, gender])), n_samples) + 
          (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))*second_mortality
      }
      
      
      transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Post TKR >=3 years < 10 years"] = (1-transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Early revision"]-
                                                                                                                             transition_matrices[i_cycle, treatment_name, , "State Post TKR <3 years", "State Death"])/4
      
      
      
      transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Post TKR >=10 years"] = (1-transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State middle revision"]-
                                                                                                                              transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Death"])/8
      
      
      # Ensure remaining patients stay in state
      # Sum probabilities of transitions to other states
      # This ensures probabilities sum to 1.
      for(i_state in 1:length(state_names)) {
        transition_matrices[i_cycle, treatment_name, , i_state, i_state] <- 1 - 
          apply(transition_matrices[i_cycle, treatment_name, , i_state, -i_state], c(1), sum, na.rm=TRUE)
      } # End loop over states
     
      
      #if(!is.null(sensitivity)) {
        #if(sensitivity == "remove_time"){
           # transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State middle revision"] <- 
           # transition_matrices[i_cycle, treatment_name, , "State Post TKR >=10 years", "State late revision"] <- 
           #transition_matrices[i_cycle, treatment_name, , "State middle revision", "State second revision"] <- 
           # transition_matrices[i_cycle, treatment_name, , "State late revision", "State second revision"] <- 
           # transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Death"]<-
           # transition_matrices[ i_cycle, treatment_name, , "State Post TKR >=10 years", "State Death"]<-
           # transition_matrices[ i_cycle, treatment_name, , "State middle revision", "State Death"]<-
           # transition_matrices[ i_cycle, treatment_name, , "State late revision", "State Death"]<-
           # transition_matrices[ i_cycle, treatment_name, , "State Post TKR <3 years", "State Post TKR >=3 years < 10 years"]<- 
           # transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Post TKR >=10 years"]<-
           # transition_matrices[i_cycle, treatment_name, , "State Post TKR >=3 years < 10 years", "State Post TKR >=3 years < 10 years"] <- 
           # transition_matrices[i_cycle, treatment_name, , "State Post TKR >=10 years", "State Post TKR >=10 years"] <- 
           # transition_matrices[i_cycle, treatment_name, , "State middle revision", "State middle revision"] <- 
            #transition_matrices[i_cycle, treatment_name, , "State late revision", "State late revision"] <- rep(0, each = n_samples)
        #}
     # }
    } # End loop over implant_names
  } # End loop over cycles
  

  i_cycle_sample <- sample(1:n_cycles, 1) 
  treatment_sample <- sample(treatment_names, 1)
  i_sample_sample <- sample(1:n_samples, 1)
  if(prod(round(rowSums(transition_matrices[i_cycle_sample,
                                            treatment_sample,
                                            i_sample_sample, , ]), digits = 10) == 1) != 1) {
    stop("Rows must sum to 1!")
  } 
  
  
  return(transition_matrices)
}
