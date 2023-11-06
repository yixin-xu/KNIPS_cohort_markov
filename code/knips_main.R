# Economic Evaluation Modelling Using R
# Advanced Markov models lecture
# Script conduct cost-effectiveness analysis of atrial fibrillation using example Markov model

# Load BCEA library to help analyse and visualise results
library(BCEA)
library(readxl)
library(ggplot2)
library(reshape)
library(dplyr)

set.seed(14142234)


source('code/generate_input_parameters.R')
source('code/generate_transition_matrices.R')
source('code/convert_transition_matrices_to_df.R')
source('code/generate_state_qalys.R')
source('code/generate_state_costs.R')
source('code/generate_net_benefit_cpp_full.R')




# Define global simulation parameters
n_samples <- 1000

# Define global model structure parameters
n_states <- 8
state_names <- paste("State", c("Post TKR <3 years", "Post TKR >=3 years < 10 years", "Post TKR >=10 years", 
                                "Early revision",  "middle revision", "late revision", "second revision", "Death"))

treatment_names <- paste("Implant", c("MoP Cem CR_Fix Mod","MoP Cem CR_Fix Mono",  "MoP Cem CR_Mob Mod",
                                      "MoP Cem PS_Fix Mod", "MoP Cem PS_Mob Mod",  "MoP Cem Con_Con Mod",
                                      "MoP Unc CR_Fix Mod", "MoP Unc CR_Mob Mod", "MoP Unc PS_Fix Mod",
                                      "MoP Hyb CR_Fix Mod", "OX Cem CR_Fix Mod", "OX Cem PS_Fix Mod"))
n_treatments <- length(treatment_names)

event_names <- c("early primary", "middle primary", "late primary", "early revision", "middle revision", "late revision", "reresion")


age_range <- "65-74"

initial_age <- substring(age_range, 1, 2)
final_age <- substring(age_range, 4, 6)


if(is.infinite(as.numeric(final_age))) {starting_age <- 90 }else{
if (initial_age == "0-") {starting_age <- 53}else{starting_age <- ceiling((as.numeric(final_age)+as.numeric(initial_age))/2)}
}
# Specify the gender
gender <- "female"

# Define global scenario parameters
if (initial_age == "0-") {n_cycles <- 50}else{
  n_cycles <- 100- as.numeric(initial_age)}

# Generate the input parameters
# This will be converted into transition matrix, state costs, and state utilities
input_parameters <- generate_input_parameters(n_samples,
                                              treatment_names = treatment_names, 
                                              state_names = state_names,
                                              initial_age = initial_age,
                                              final_age = final_age,
                                              starting_age = starting_age,
                                              gender = gender)

# Run the Markov model to get the model outputs
model_outputs <- generate_net_benefit(input_parameters,
                                      treatment_names = treatment_names, 
                                      state_names = state_names,
                                      initial_age = initial_age,
                                      final_age = final_age,
                                      starting_age = starting_age,
                                      gender = gender,
                                      n_cycles = n_cycles,
                                      lambda = 20000)

####################################################################################################
## Quick manual check of resutls ###################################################################
####################################################################################################

# Expected net benefit at ?20,000
#with(model_outputs, colMeans(20000 * total_qalys - total_costs))

#rowMeans(model_outputs$incremental_net_benefit)
##################################################################################################################
## Now use BCEA to analyse the outputs   #########################################################################
##################################################################################################################
# Create a bcea object for the knips model
# Note costs and QALYs need to be transposed in this example for BCEA to run
knips_bcea <- bcea(e = t(model_outputs$total_qalys), c = t(model_outputs$total_costs), ref = 1, interventions = treatment_names) 
# Get summary statistics
#summary(knips_bcea, wtp = 20000)

# Plot a CEAC
setwd('C:/Users/yx18392/Desktop/knips/cohort markov model/results')
png(file=paste0(gender,"-", age_range,"_multice.png"))
knips_multi_ce <- multi.ce(knips_bcea)
ceac.plot(knips_multi_ce, graph = "ggplot",
         line = list(colors = rainbow(12)),
         pos = c(1,1))


dev.off()

format_results <- function(x, n_digits = 2) {
  paste0(format(mean(x), nsmall = n_digits, digits = n_digits), " (",
         format(quantile(x, prob = 0.025), nsmall = n_digits, digits = n_digits), " ,",
         format(quantile(x, prob = 0.975), nsmall = n_digits, digits = n_digits), ")")
}

results_matrix <- matrix(nrow = n_treatments, ncol = 4)
rownames(results_matrix) <- treatment_names
for(treatment in treatment_names) {
  results_matrix[treatment, ] <- c(format_results(model_outputs$total_qalys[treatment,]),
                                   format_results(model_outputs$total_costs[treatment,]),
                                   format_results(model_outputs$net_benefit[treatment,]),
                                   format_results(model_outputs$incremental_net_benefit[treatment,]))
}

results_matrix = cbind(treatment_names,results_matrix, round(rowMeans(model_outputs$ICER),2))
colnames(results_matrix)  = c("Treatment_names", "QALYs", "Costs", "Net_benefit", "Inc.Net_benefit", "ICER")

for (treatment in treatment_names){
  if(results_matrix[treatment,6]<0) {
    if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[1,])>0){results_matrix[treatment,6] = "Dominant"}
    if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[1,])<0){results_matrix[treatment,6] = "Dominated"} 
  }
  
}

icer_table = data.frame(results_matrix)
write.csv(icer_table, file = paste0(gender,"-", age_range,"icer_results.csv"))


##################################################################################
## VoI analysis - takes a few minutes ############################################
##################################################################################

evpi_table <- matrix(nrow = 7, ncol = 2)
rownames(evpi_table) <- c("Total", "1st revision probabilities_<3", "1st revision probabilities_3-10",
                          "1st revision probabilities_>10","2nd and higher revision probabilities",
                          "primary utilities", "revision utilities")
colnames(evpi_table) <- c("Per person", "Population")

# Assume implant decision remains relevant for at least 10 years
# And that there are 160000 THR per year
# https://www.njrcentre.org.uk/njrcentre/Patients/Joint-replacement-statistics#:~:text=In%20England%20and%20Wales%20there,the%20practice%20is%20growing%20rapidly.
# This is regardless of age and gender but assume EVPI same for each category
technology_horizon <- 10
discounted_population_size <- sum((1/1.035)^(0:(technology_horizon - 1))) * 0.2218 * 92185 

# Total EVPI
evpi_table["Total", c("Per person", "Population")] <-  knips_bcea$evi[201] * c(1, discounted_population_size)


# Log rate of first revision (same EVPPI as if calculating probability)
# Use GP instead of GAM as 4 parameters
evppi_gp_1st_revision_1 <- evppi(param_idx = c("log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mono",
                                   "log_rate_1st_revision_<3Implant MoP Cem CR_Fix Mod",
                                   "log_rate_1st_revision_<3Implant MoP Cem CR_Mob Mod",
                                   "log_rate_1st_revision_<3Implant MoP Cem PS_Fix Mod",
                                   "log_rate_1st_revision_<3Implant MoP Cem PS_Mob Mod",
                                   "log_rate_1st_revision_<3Implant MoP Cem Con_Con Mod",
                                   "log_rate_1st_revision_<3Implant MoP Unc CR_Fix Mod",
                                   "log_rate_1st_revision_<3Implant MoP Unc CR_Mob Mod",
                                   "log_rate_1st_revision_<3Implant MoP Unc PS_Fix Mod",
                                   "log_rate_1st_revision_<3Implant MoP Hyb CR_Fix Mod",
                                   "log_rate_1st_revision_<3Implant OX Cem CR_Fix Mod",
                                   "log_rate_1st_revision_<3Implant OX Cem PS_Fix Mod"),
                                 input = input_parameters, he = knips_bcea, method = 'gp')
evpi_table["1st revision probabilities_<3", c("Per person", "Population")] <- evppi_gp_1st_revision_1$evppi[201] * c(1, discounted_population_size)

evppi_gp_1st_revision_2 <- evppi(param_idx = c("log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mono",
                                                "log_rate_1st_revision_3-10Implant MoP Cem CR_Fix Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Cem CR_Mob Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Cem PS_Fix Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Cem PS_Mob Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Cem Con_Con Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Unc CR_Fix Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Unc CR_Mob Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Unc PS_Fix Mod",
                                                "log_rate_1st_revision_3-10Implant MoP Hyb CR_Fix Mod",
                                                "log_rate_1st_revision_3-10Implant OX Cem CR_Fix Mod",
                                                "log_rate_1st_revision_3-10Implant OX Cem PS_Fix Mod"),
                                  input = input_parameters, he = knips_bcea, method = 'gp')
evpi_table["1st revision probabilities_3-10", c("Per person", "Population")] <- evppi_gp_1st_revision_2$evppi[201] * c(1, discounted_population_size)

evppi_gp_1st_revision_3 <- evppi(param_idx = c("log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mono",
                                               "log_rate_1st_revision_>10Implant MoP Cem CR_Fix Mod",
                                               "log_rate_1st_revision_>10Implant MoP Cem CR_Mob Mod",
                                               "log_rate_1st_revision_>10Implant MoP Cem PS_Fix Mod",
                                               "log_rate_1st_revision_>10Implant MoP Cem PS_Mob Mod",
                                               "log_rate_1st_revision_>10Implant MoP Cem Con_Con Mod",
                                               "log_rate_1st_revision_>10Implant MoP Unc CR_Fix Mod",
                                               "log_rate_1st_revision_>10Implant MoP Unc CR_Mob Mod",
                                               "log_rate_1st_revision_>10Implant MoP Unc PS_Fix Mod",
                                               "log_rate_1st_revision_>10Implant MoP Hyb CR_Fix Mod",
                                               "log_rate_1st_revision_>10Implant OX Cem CR_Fix Mod",
                                               "log_rate_1st_revision_>10Implant OX Cem PS_Fix Mod"),
                                 input = input_parameters, he = knips_bcea, method = 'gp')
evpi_table["1st revision probabilities_>10", c("Per person", "Population")] <- evppi_gp_1st_revision_3$evppi[201] * c(1, discounted_population_size)


# Log rate 2nd and higher revision probabilities
evppi_gam_2nd_and_higher_revision <- evppi(param_idx =  
                                             c("log_rate_2nd_revision_early", "log_rate_2nd_revision_middle", 
                                               "log_rate_2nd_revision_late", "log_rate_higher_revision"), 
                                           input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["2nd and higher revision probabilities", c("Per person", "Population")] <- evppi_gam_2nd_and_higher_revision$evppi[201] * c(1, discounted_population_size)

# Utilities
evppi_gam_primary_utilities <- evppi(param_idx =  
                               c( "qalys_State Post TKR <3 years", 
                                  "qalys_State Post TKR >=3 years < 10 years", 
                                  "qalys_State Post TKR >=10 years"), 
                             input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["primary utilities", c("Per person", "Population")] <- evppi_gam_primary_utilities$evppi[201] * c(1, discounted_population_size)

evppi_gam_revision_utilities <- evppi(param_idx =  
                               c( "qalys_State Early revision", 
                                  "qalys_State middle revision", 
                                  "qalys_State late revision", 
                                  "qalys_State second revision"), 
                             input = input_parameters, he = knips_bcea, method = 'gam')
evpi_table["revision utilities", c("Per person", "Population")] <- evppi_gam_revision_utilities$evppi[201] * c(1, discounted_population_size)




