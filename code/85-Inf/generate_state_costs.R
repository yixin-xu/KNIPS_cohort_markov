generate_state_costs<-function(input_parameters,
                               treatment_names = treatment_names, 
                               state_names = state_names,
                               initial_age = initial_age,
                               final_age = final_age,
                               starting_age = starting_age,
                               gender = gender,
                               sensitivity=NULL) {
  
  n_treatments <- length(treatment_names)
  n_samples <- dim(input_parameters)[1]
  n_states <- length(state_names)
  
  # Construct array of state costs with default value of zero
  state_costs <- array(0, dim = c(n_samples, n_treatments, n_states),
                       dimnames = list( NULL, treatment_names, state_names))
  
  # State costs are a mixture of treatment costs, acute event costs, event managment costs, 
  # and costs of transient events
  
  for(treatment_name in treatment_names) {
    state_costs[ ,treatment_name , "State Post TKR <3 years"] <- input_parameters[ , "cost_revision"] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))
    state_costs[ ,treatment_name , "State Post TKR >=3 years < 10 years"] <- input_parameters[ , "cost_revision"] * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))
    
    
  }
  
  state_costs[ , , "State Early revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))
  state_costs[ , , "State middle revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))
 
  state_costs[ , , "State second revision"] <- input_parameters[ , "cost_revision"] * (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))
  
  if(!is.null(sensitivity)) {
    if(sensitivity == "remove_2nd_higher_rate") {
      state_costs[ , , "State second revision"] = rep(rep(0,n = n_samples), n= n_treatments)}
    if(sensitivity == "remove_time") { 
      state_costs[ , , "State Post TKR >=3 years < 10 years"] <-state_costs[ , , "State Post TKR >=10 years"] <-
        state_costs[ , , "State middle revision"] <-  state_costs[ , , "State late revision"] <-rep(rep(0,n = n_samples), n= n_treatments)
    }
  }
  
  # End loop over treatments
  
  return(state_costs)
} # End function
