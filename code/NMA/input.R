generate_input_parameters <- function(n_samples, treatment_names = treatment_names, 
                                      state_names = state_names,
                                      initial_age = initial_age,
                                      final_age = final_age,
                                      starting_age = starting_age,
                                      gender = gender,
                                      sensitivity = NULL) {
  
  n_treatments <- length(treatment_names)
  
  
  parameter_names <- c(paste0("log_rate_1st_revision_<3", treatment_names),
                       paste0("log_rate_1st_revision_3-10", treatment_names),
                       paste0("log_rate_1st_revision_>10", treatment_names),
                       "log_rate_2nd_revision_early", "log_rate_2nd_revision_middle", 
                       "log_rate_2nd_revision_late", "log_rate_higher_revision",
                       "cost_revision",
                       paste0("qalys_", state_names),
                       paste0("implant_cost_", treatment_names))
  n_parameters <- length(parameter_names)
  input_parameters <- array(dim = c(n_samples, n_parameters), dimnames = list(NULL, parameter_names))
  
  # data: female 55-64 years old group
  lifetables <- read_excel(paste0(data_directory, "/KNIPS Main input data.xlsx"), sheet = "uk_lifetables")
  lograte_revision <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "revision_log_rate")
  costs <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "costs")
  utilities <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities")
  un_utilities <- read_excel(paste0(data_directory, "/cohort_model_inputs.xlsx"), sheet = "utilities_unadjusted")
  log_rate_1st_revision <- read_excel(paste0(data_directory,"/", paste0(gender, "-", initial_age, "-", "rate.xlsx")),sheet = "Sheet4" )
  
  
  # Impute 1st and 12th implants to be the average over all implants
  
  if ("log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[1], sd = (log_rate_1st_revision$UL1[1]-log_rate_1st_revision$LL1[1])/2*1.96)}
  
  
  if ("log_rate_1st_revision_<3Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem CR_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[2], sd = (log_rate_1st_revision$UL1[2]-log_rate_1st_revision$LL1[2])/2*1.96)}                                                                                      
  
  if ("log_rate_1st_revision_<3Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem PS_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[3], sd = (log_rate_1st_revision$UL1[3]-log_rate_1st_revision$LL1[3])/2*1.96)}
  
  if ("log_rate_1st_revision_<3Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Cem PS_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[4], sd = (log_rate_1st_revision$UL1[4]-log_rate_1st_revision$LL1[4])/2*1.96)}                                                                                      
  
  if ("log_rate_1st_revision_<3Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Unc CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[5], sd = (log_rate_1st_revision$UL1[5]-log_rate_1st_revision$LL1[5])/2*1.96)}                                                                                       
  
  if ("log_rate_1st_revision_<3Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Unc CR_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[6], sd = (log_rate_1st_revision$UL1[6]-log_rate_1st_revision$LL1[6])/2*1.96)}                                                                                      
  
  if ("log_rate_1st_revision_<3Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_<3Implant MoP Hyb CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean1[7], sd = (log_rate_1st_revision$UL1[7]-log_rate_1st_revision$LL1[7])/2*1.96)}                                                                                       
  
  if ("log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[1], sd = (log_rate_1st_revision$UL2[1]-log_rate_1st_revision$LL2[1])/2*1.96)}                                                                                       
  if ("log_rate_1st_revision_3-10Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem CR_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[2], sd = (log_rate_1st_revision$UL2[2]-log_rate_1st_revision$LL2[2])/2*1.96)}                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem PS_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[3], sd = (log_rate_1st_revision$UL2[3]-log_rate_1st_revision$LL2[3])/2*1.96)}                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Cem PS_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[4], sd = (log_rate_1st_revision$UL2[4]-log_rate_1st_revision$LL2[4])/2*1.96)}                                                                                       
  if ("log_rate_1st_revision_3-10Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Unc CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[5], sd = (log_rate_1st_revision$UL2[5]-log_rate_1st_revision$LL2[5])/2*1.96)}                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Unc CR_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[6], sd = (log_rate_1st_revision$UL2[6]-log_rate_1st_revision$LL2[6])/2*1.96)}                                                                                      
  if ("log_rate_1st_revision_3-10Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_3-10Implant MoP Hyb CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean2[7], sd = (log_rate_1st_revision$UL2[7]-log_rate_1st_revision$LL2[7])/2*1.96)}                                                                                      
  
  
  
  if ("log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[1], sd = (log_rate_1st_revision$UL3[1]-log_rate_1st_revision$LL3[1])/2*1.96)}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem CR_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[2], sd = (log_rate_1st_revision$UL3[2]-log_rate_1st_revision$LL3[2])/2*1.96)}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem PS_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[3], sd = (log_rate_1st_revision$UL3[3]-log_rate_1st_revision$LL3[3])/2*1.96)}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Cem PS_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[4], sd = (log_rate_1st_revision$UL3[4]-log_rate_1st_revision$LL3[4])/2*1.96)}                                                                                      
  if ("log_rate_1st_revision_>10Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Unc CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[5], sd = (log_rate_1st_revision$UL3[5]-log_rate_1st_revision$LL3[5])/2*1.96)}                                                                                       
  if ("log_rate_1st_revision_>10Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Unc CR_Mob Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[6], sd = (log_rate_1st_revision$UL3[6]-log_rate_1st_revision$LL3[6])/2*1.96) }                                                                                      
  if ("log_rate_1st_revision_>10Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "log_rate_1st_revision_>10Implant MoP Hyb CR_Fix Mod"] <- rnorm(n_samples, mean = log_rate_1st_revision$mean3[7], sd = (log_rate_1st_revision$UL3[7]-log_rate_1st_revision$LL3[7])/2*1.96)}                                                                                       
  
  lograte_revision[,2:4]= log(lograte_revision[,2:4])
  
  input_parameters[ , "log_rate_2nd_revision_early"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "early_revision"],
                                                                                     sd = ((UL[parameter == "early_revision"]-LL[parameter == "early_revision"])/2*1.96)))
  input_parameters[ , "log_rate_2nd_revision_middle"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "middle_revision"],
                                                                                      sd = ((UL[parameter == "middle_revision"]-LL[parameter == "middle_revision"])/2*1.96)))
  
  if(is.infinite(as.numeric(final_age))) {
    input_parameters[ , "log_rate_2nd_revision_late"] <- NaN}else{
      input_parameters[ , "log_rate_2nd_revision_late"] <- with(lograte_revision, rnorm(n_samples, mean =  mean[parameter == "late_revision"],
                                                                                        sd = ((UL[parameter == "late_revision"]-LL[parameter == "late_revision"])/2*1.96)))
      
    }
  
  
  # third revision
  sd_third = (as.numeric(lograte_revision[5,4]) - as.numeric(lograte_revision[5,3]))/2*1.96
  se_third = sd_third/sqrt(as.numeric(lograte_revision[5,5]))
  tau_third = 1/se_third^2
  
  # fourth revision
  sd_fourth = (as.numeric(lograte_revision[6,4]) - as.numeric(lograte_revision[6,3]))/2*1.96
  se_fourth = sd_fourth/sqrt(as.numeric(lograte_revision[6,5]))
  tau_fourth = 1/se_fourth^2
  
  # fifth revision
  sd_fifth = (as.numeric(lograte_revision[7,4]) - as.numeric(lograte_revision[7,3]))/2*1.96
  se_fifth = sd_fifth/sqrt(as.numeric(lograte_revision[7,5]))
  tau_fifth = 1/se_fifth^2
  
  # sixth revision
  sd_sixth = (as.numeric(lograte_revision[8,4]) - as.numeric(lograte_revision[8,3]))/2*1.96
  se_sixth = sd_sixth/sqrt(as.numeric(lograte_revision[8,5]))
  tau_sixth = 1/se_sixth^2
  
  # seventh revision
  sd_seventh = (as.numeric(lograte_revision[9,4]) - as.numeric(lograte_revision[9,3]))/2*1.96
  se_seventh = sd_seventh/sqrt(as.numeric(lograte_revision[9,5]))
  tau_seventh = 1/se_seventh^2
  
  
  # weight of ecah revision, using in the calculation of weighted average of log rates
  w_third = tau_third/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh )
  w_fourth = tau_fourth/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  w_fifth = tau_fifth/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  w_sixth = tau_sixth/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  w_seventh = tau_seventh/(tau_third + tau_fourth + tau_fifth + tau_sixth + tau_seventh)
  
  average_lograte = w_third*as.numeric(lograte_revision[5,3]) + w_fourth*as.numeric(lograte_revision[6,3]) + w_fifth*as.numeric(lograte_revision[7,3]) + 
    w_sixth*as.numeric(lograte_revision[8,3]) + w_seventh*as.numeric(lograte_revision[9,3])
  average_se = sqrt(w_third^2*se_third^2 + w_fourth^2*se_fourth^2 + w_fifth^2*se_fifth^2 + w_sixth^2*se_sixth^2 + w_seventh^2*se_seventh^2)
  
  input_parameters[ , "log_rate_higher_revision"] = rnorm(n_samples, average_lograte, average_se)
  
  
  if ("implant_cost_Implant MoP Cem CR_Fix Mono" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem CR_Fix Mono"] <- rnorm(n_samples, as.numeric(costs[1,2]), sd = ((as.numeric(costs[1,5])-as.numeric(costs[1,4]))/2*1.96)) }
  if ("implant_cost_Implant MoP Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem CR_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[2,2]), sd = ((as.numeric(costs[2,5])-as.numeric(costs[2,4]))/2*1.96)) }
  if ("implant_cost_Implant MoP Cem CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem CR_Mob Mod"] <-  rnorm(n_samples, as.numeric(costs[3,2]), sd = ((as.numeric(costs[3,5])-as.numeric(costs[3,4]))/2*1.96)) }
  if ("implant_cost_Implant MoP Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem PS_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[4,2]), sd = ((as.numeric(costs[4,5])-as.numeric(costs[4,4]))/2*1.96))}
  if ("implant_cost_Implant MoP Cem PS_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem PS_Mob Mod"] <- rnorm(n_samples, as.numeric(costs[5,2]), sd = ((as.numeric(costs[5,5])-as.numeric(costs[5,4]))/2*1.96))}
  if ("implant_cost_Implant MoP Cem Con_Con Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Cem Con_Con Mod"] <- rnorm(n_samples, as.numeric(costs[6,2]), sd = ((as.numeric(costs[6,5])-as.numeric(costs[6,4]))/2*1.96)) }
  if ("implant_cost_Implant MoP Unc CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Unc CR_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[7,2]), sd = ((as.numeric(costs[7,5])-as.numeric(costs[7,4]))/2*1.96))}
  if ("implant_cost_Implant MoP Unc CR_Mob Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Unc CR_Mob Mod"] <- rnorm(n_samples, as.numeric(costs[8,2]), sd = ((as.numeric(costs[8,5])-as.numeric(costs[8,4]))/2*1.96)) }
  if ("implant_cost_Implant MoP Unc PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Unc PS_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[9,2]), sd = ((as.numeric(costs[9,5])-as.numeric(costs[9,4]))/2*1.96)) }
  if ("implant_cost_Implant MoP Hyb CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant MoP Hyb CR_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[10,2]), sd = ((as.numeric(costs[10,5])-as.numeric(costs[10,4]))/2*1.96)) }
  if ("implant_cost_Implant OX Cem CR_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant OX Cem CR_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[11,2]), sd = ((as.numeric(costs[11,5])-as.numeric(costs[11,4]))/2*1.96)) }
  if ("implant_cost_Implant OX Cem PS_Fix Mod" %in% colnames(input_parameters)){
    input_parameters[ , "implant_cost_Implant OX Cem PS_Fix Mod"] <- rnorm(n_samples, as.numeric(costs[12,2]), sd = ((as.numeric(costs[12,5])-as.numeric(costs[12,4]))/2*1.96))  }
  
  
  input_parameters[ , "cost_revision"] <- rep(as.numeric(costs[15,2]), n = n_samples) 
  
  # states qalys
  if (ini_age == "50") {ini_age2 <- 0}else{
    ini_age2 <- ini_age}
  utility_row_index <- which(utilities[, "age"] == ini_age2 &
                               utilities[, "gender"] == gender)
  
  
  input_parameters[ ,"qalys_State Post TKR <3 years"] <- input_parameters[ ,"qalys_State Post TKR >=3 years < 10 years"] <- input_parameters[ ,"qalys_State Post TKR >=10 years"] <-
    rep(rnorm(n_samples, mean = as.numeric(utilities[utility_row_index,6]), sd = (as.numeric(utilities[utility_row_index,8])-as.numeric(utilities[utility_row_index,7]))/2*1.96), n= n_treatments)
  input_parameters[ ,"qalys_State Early revision"] <- input_parameters[ ,"qalys_State middle revision"] <-  input_parameters[ ,"qalys_State late revision"] <- input_parameters[ ,"qalys_State second revision"] <-
    rep(rnorm(n_samples, mean = as.numeric(utilities[utility_row_index,12]), sd = (as.numeric(utilities[utility_row_index,14])-as.numeric(utilities[utility_row_index,13]))/2*1.96), n= n_treatments)
  input_parameters[ ,"qalys_State Death"]<- rep(rep(0, each = n_samples), n= n_treatments)
  
  
  #
  if(!is.null(sensitivity)) {
    if(sensitivity == "un_utilities") {
      input_parameters[ ,"qalys_State Post TKR <3 years"] <- input_parameters[ ,"qalys_State Post TKR >=3 years < 10 years"] <- input_parameters[ ,"qalys_State Post TKR >=10 years"] <- 
        rep(rnorm(n_samples, mean = as.numeric(un_utilities[utility_row_index,6]), sd = (as.numeric(un_utilities[utility_row_index,8])-as.numeric(un_utilities[utility_row_index,7]))/2*1.96), n= n_treatments)
      input_parameters[ ,"qalys_State Early revision"] <- input_parameters[ ,"qalys_State middle revision"] <-  input_parameters[ ,"qalys_State late revision"] <- input_parameters[ ,"qalys_State second revision"] <-
        rep(rnorm(n_samples, mean = as.numeric(un_utilities[utility_row_index,12]), sd = (as.numeric(un_utilities[utility_row_index,14])-as.numeric(un_utilities[utility_row_index,13]))/2*1.96), n= n_treatments)
      input_parameters[ ,"qalys_State Death"]<- rep(rep(0, each = n_samples), n= n_treatments)
    }
    if(sensitivity == "remove_2nd_higher_rate"){
      input_parameters[ ,"log_rate_higher_revision"] <- rep(rep(NaN, each = n_samples), n= n_treatments)
      input_parameters[ ,"qalys_State second revision"] <-rep(rep(0, each = n_samples), n= n_treatments)
    }
    if(sensitivity == "remove_time"){
      input_parameters[, grepl("log_rate_1st_revision_3-10", colnames(input_parameters))] = input_parameters[, grepl("log_rate_1st_revision_>10", colnames(input_parameters))] = 
        input_parameters[, "log_rate_2nd_revision_middle"] = input_parameters[, "log_rate_2nd_revision_late"] = rep(rep(NaN, each = n_samples), n= n_treatments)
      input_parameters[ ,"qalys_State Post TKR >=3 years < 10 years"] <- input_parameters[ ,"qalys_State Post TKR >=10 years"] <-
        input_parameters[ ,"qalys_State middle revision"] <-  input_parameters[ ,"qalys_State late revision"] <- rep(rep(0, each = n_samples), n= n_treatments)
    }
  }
  
  
  
  return(input_parameters)
} # End function

