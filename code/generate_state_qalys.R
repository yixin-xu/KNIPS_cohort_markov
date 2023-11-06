# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Function to convert input parameters state utilities

generate_state_qalys <- function(input_parameters, 
                                 treatment_names = treatment_names, 
                                 state_names = state_names,
                                 initial_age = initial_age,
                                 final_age = final_age,
                                 starting_age = starting_age,
                                 gender = gender,
                                 sensitivity = NULL) {
  
  n_treatments <- length(treatment_names)
  n_samples <- dim(input_parameters)[1]
  n_states <- length(state_names)
  
  utilities <- read_excel(paste0(data_directory, "/cohort model inputs.xlsx"), sheet = "utilities")
  un_utilities <- read_excel(paste0(data_directory, "/cohort model inputs.xlsx"), sheet = "utilities_unadjusted")
 
  # Construct array of state utilities with default value 0
  state_utilities <- array(dim = c(n_samples, n_treatments, n_states), dimnames = list(NULL, treatment_names,state_names))
  event_disutilities <- array(dim = c(n_samples, n_treatments, n_states), dimnames = list(NULL, treatment_names,state_names))
  # State utilities are a mixture of state event  utilities, and disutilities of transient events
  
  
  multi_disutilities = rlnorm(n_samples, log(4.5), log(1.5))
  
  norm_disutilities <- rnorm(n_samples, mean = as.numeric(utilities[which(grepl(paste0(initial_age," ", gender),utilities$...1)),"disutilities"]),  
                             sqrt((as.numeric(utilities[which(grepl(paste0(initial_age," ", gender),utilities$...1)),"SE_pre"])^2+(as.numeric(utilities[which(grepl(paste0(initial_age," ", gender),utilities$...1)),"SE_ad_6 months after"])^2))))
  
  
  
  if(!is.null(sensitivity)) {
    if(sensitivity == "un_utilities") {
      norm_disutilities <- rnorm(n_samples, mean = as.numeric(un_utilities[which(grepl(paste0(initial_age," ", gender),un_utilities$...1)),14]),  sqrt((as.numeric(un_utilities[which(grepl(paste0(initial_age," ", gender),un_utilities$...1)),17])^2+(as.numeric(un_utilities[which(grepl(paste0(initial_age," ", gender),un_utilities$...1)),18])^2))))
    }
  }
  for(treatment_name in treatment_names) {
    event_disutilities[ ,treatment_name , "State Post TKR <3 years"] <- norm_disutilities * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_<3", treatment_name)])))* multi_disutilities/2
    event_disutilities[ ,treatment_name , "State Post TKR >=3 years < 10 years"] <- norm_disutilities * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_3-10", treatment_name)])))* multi_disutilities/2
    event_disutilities[ ,treatment_name , "State Post TKR >=10 years"] <- norm_disutilities * 
      (1 - exp(-exp(input_parameters[, paste0("log_rate_1st_revision_>10", treatment_name)])))* multi_disutilities/2
    
  }
  
  
  event_disutilities[ , , "State Early revision"] <- norm_disutilities * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_early"])))* multi_disutilities/2
  event_disutilities[ , , "State middle revision"] <- norm_disutilities * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_middle"])))* multi_disutilities/2
  event_disutilities[ , , "State late revision"] <- norm_disutilities * (1 - exp(-exp(input_parameters[ , "log_rate_2nd_revision_late"])))* multi_disutilities/2
  event_disutilities[ , , "State second revision"] <- norm_disutilities * (1 - exp(-exp(input_parameters[ , "log_rate_higher_revision"])))* multi_disutilities/2
  
  
  state_utilities[, , "State Post TKR <3 years"] = rep(input_parameters[ ,"qalys_State Post TKR <3 years"], n= n_treatments) + event_disutilities[, , "State Post TKR <3 years"]
  state_utilities[, , "State Post TKR >=3 years < 10 years"] = rep(input_parameters[ ,"qalys_State Post TKR >=3 years < 10 years"], n= n_treatments) + event_disutilities[, , "State Post TKR >=3 years < 10 years"]
  state_utilities[, , "State Post TKR >=10 years"] =rep(input_parameters[ ,"qalys_State Post TKR >=10 years"], n= n_treatments) + event_disutilities[, , "State Post TKR >=10 years"]
  state_utilities[, , "State Early revision"] = rep(input_parameters[ ,"qalys_State Early revision"], n= n_treatments) + event_disutilities[, , "State Early revision"]
  state_utilities[, , "State middle revision"] = rep(input_parameters[ ,"qalys_State middle revision"], n= n_treatments) + event_disutilities[, , "State middle revision"]
  state_utilities[, , "State late revision"] = rep(input_parameters[ ,"qalys_State late revision"], n= n_treatments) + event_disutilities[, , "State late revision"]
  state_utilities[, , "State second revision"] = rep(input_parameters[ ,"qalys_State second revision"], n= n_treatments) + event_disutilities[, , "State second revision"]
  state_utilities[, , "State Death"] = rep(input_parameters[ ,"qalys_State Death"], n= n_treatments)
  
  if(!is.null(sensitivity)) {
    if(sensitivity == "remove_2nd_higher_rate") {
      state_utilities[, , "State second revision"] = rep(input_parameters[ ,"qalys_State second revision"], n= n_treatments)}
    if(sensitivity == "remove_time"){ 
      state_utilities[, , "State Post TKR >=3 years < 10 years"] = state_utilities[, , "State Post TKR >=10 years"] =    
        state_utilities[, , "State middle revision"] =  state_utilities[, , "State late revision"] =  rep(rep(0,n = n_samples), n= n_treatments)
      
    }
  }
 
  return(state_utilities)
} # End function
