library(BCEA)
library(readxl)
library(ggplot2)
library(reshape)
library(dplyr)

set.seed(14142234)


source('code/generate_input_parameters.R')
source('code/generate_transition_matrices.R')
source('code/generate_state_qalys.R')
source('code/generate_state_costs.R')
source('code/generate_net_benefit_lapply.R')


data_directory = "C:/Users/yx18392/Desktop/semi Markov data"

# Define global simulation parameters
n_samples <- 1000

# Define global model structure parameters
n_states <- 8
state_names <- paste("State", c("Post TKR <3 years", "Post TKR >=3 years < 10 years", "Post TKR >=10 years", 
                                "Early revision",  "middle revision", "late revision", "second revision", "Death"))

treatment_names <- paste("Implant", c("MoP Cem CR_Fix Mod","MoP Cem CR_Mob Mod",  "MoP Cem PS_Fix Mod",
                                      "MoP Cem PS_Mob Mod", "MoP Unc CR_Fix Mod",  "MoP Unc CR_Mob Mod",
                                      "MoP Hyb CR_Fix Mod"))
n_treatments <- length(treatment_names)

event_names <- c("early primary", "middle primary", "late primary", "early revision", "middle revision", "late revision", "reresion")


age_range <- "85-Inf"

initial_age <- substring(age_range, 1, 2)
final_age <- substring(age_range, 4, 6)
if (final_age == "5") {final_age <- 55}

if (initial_age == "0-") {ini_age <- 50}else{
  ini_age <- as.numeric(initial_age)}


if(is.infinite(as.numeric(final_age))) {starting_age <- 90 }else{
  if (initial_age == "0-") {starting_age <- 53}else{starting_age <- ceiling((as.numeric(final_age)+as.numeric(initial_age))/2)}
}
# Specify the gender
gender <- "male"

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
model_outputs <- generate_net_benefit_lapply(input_parameters, lambda = 20000)

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
setwd('C:/Users/yx18392/Desktop')

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

results_matrix = cbind(treatment_names,results_matrix, round(model_outputs$ICER,2))
colnames(results_matrix)  = c("Treatment_names", "QALYs", "Costs", "Net_benefit", "Inc.Net_benefit", "ICER")

for (treatment in treatment_names){
  if(results_matrix[treatment,6]<0) {
    if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[1,])>0){results_matrix[treatment,6] = "Dominant"}
    if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[1,])<0){results_matrix[treatment,6] = "Dominated"} 
  }
  
}

icer_table = data.frame(results_matrix)
write.csv(icer_table, file = paste0(gender,"-", age_range,"icer_results.csv"))


