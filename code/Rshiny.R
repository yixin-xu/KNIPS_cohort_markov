library(BCEA)
library(readxl)
library(ggplot2)
library(reshape)
library(dplyr)
library(shiny)

set.seed(14142234)


source('code/generate_input_parameters.R')
source('code/generate_transition_matrices.R')
source('code/convert_transition_matrices_to_df.R')
source('code/generate_state_qalys.R')
source('code/generate_state_costs.R')
source('code/generate_net_benefit_cpp_full.R')

treatment_names_global <- list("Implant MoP Cem CR_Fix Mod","Implant MoP Cem CR_Fix Mono",  "Implant MoP Cem CR_Mob Mod",
                               "Implant MoP Cem PS_Fix Mod", "Implant MoP Cem PS_Mob Mod",  "Implant MoP Cem Con_Con Mod",
                               "Implant MoP Unc CR_Fix Mod", "Implant MoP Unc CR_Mob Mod", "Implant MoP Unc PS_Fix Mod",
                               "Implant MoP Hyb CR_Fix Mod", "Implant OX Cem CR_Fix Mod", "Implant OX Cem PS_Fix Mod")
n_treatments_global <- length(treatment_names_global)

#n_states <- 8
#state_names <- paste("State", c("Post TKR <3 years", "Post TKR >=3 years < 10 years", "Post TKR >=10 years", 
#                                "Early revision",  "middle revision", "late revision", "second revision", "Death"))





f_wrapper <- function(age_range, gender, treatment_names, n_samples) {
  
  treatment_names = unlist(treatment_names)
  gender = as.character(gender)
  age_range = as.character(age_range)
  
  initial_age = substring(age_range, 1, 2)
  final_age = substring(age_range, 4, 6)
  
  
  if(is.infinite(as.numeric(final_age))) {starting_age <- 90 }else{
    if (initial_age == "0-") {starting_age <- 53}else{starting_age <- ceiling((as.numeric(final_age)+as.numeric(initial_age))/2)}
  }
  
  if (initial_age == "0-") {n_cycles <- 50}else{
    n_cycles <- 100- as.numeric(initial_age)}
  
  n_treatments <- length(treatment_names)
  
  n_states <- 8
  state_names <- paste("State", c("Post TKR <3 years", "Post TKR >=3 years < 10 years", "Post TKR >=10 years", 
                                  "Early revision",  "middle revision", "late revision", "second revision", "Death"))
  
  
  
  
  input_parameters <- generate_input_parameters(n_samples, treatment_names = treatment_names, 
                                                state_names = state_names,
                                                initial_age = initial_age,
                                                final_age = final_age,
                                                starting_age = starting_age,
                                                gender = gender,
                                                sensitivity = NULL)
  
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
  
  return(model_outputs) 
}


ui  <- fluidPage (#creates empty page
  
  # title of app
  titlePanel("Knips cohort Markov Model in Shiny"),
  #layout is a sidebar—layout
  sidebarLayout(
    sidebarPanel(#open sidebar panel
      #input type 
      selectInput(inputId ="SI_gender",label ="Gender",
                  choices = list("female", "male"), 
                  selected = "female"),
      selectInput(inputId ="SI_age_range",label ="Age range",
                  choices = list("0-55", "55-64","65-74", "75-84", "85-Inf"),
                  selected = "55-64"),
      checkboxGroupInput(inputId ="SI_treatment_names",
                         label = "Treatment name",
                         choices = list("Implant MoP Cem CR_Fix Mod","Implant MoP Cem CR_Fix Mono",  "Implant MoP Cem CR_Mob Mod",
                                        "Implant MoP Cem PS_Fix Mod", "Implant MoP Cem PS_Mob Mod",  "Implant MoP Cem Con_Con Mod",
                                        "Implant MoP Unc CR_Fix Mod", "Implant MoP Unc CR_Mob Mod", "Implant MoP Unc PS_Fix Mod",
                                        "Implant MoP Hyb CR_Fix Mod", "Implant OX Cem CR_Fix Mod", "Implant OX Cem PS_Fix Mod"),
                         selected = list("Implant MoP Cem CR_Fix Mod","Implant MoP Cem CR_Fix Mono",  "Implant MoP Cem CR_Mob Mod",
                                         "Implant MoP Cem PS_Fix Mod", "Implant MoP Cem PS_Mob Mod",  "Implant MoP Cem Con_Con Mod",
                                         "Implant MoP Unc CR_Fix Mod", "Implant MoP Unc CR_Mob Mod", "Implant MoP Unc PS_Fix Mod",
                                         "Implant MoP Hyb CR_Fix Mod", "Implant OX Cem CR_Fix Mod", "Implant OX Cem PS_Fix Mod"),
                         inline = TRUE),
      numericInput(inputId ="SI_n_samples",
                   label ="Number of samples",
                   value = 100,
                   min= 0,
                   max= 1000),
      #action button runs model when pressed
      actionButton(inputId ="run_model",
                   label   ="Run model")
      
    ),#close sidebarPanel
    #open main panel
    mainPanel(
      tabsetPanel(
       tabPanel(
          title = "Results Table",
          tableOutput(outputId ="SO_icer_table")
        ),
        tabPanel(
          title = "Cost-effectiveness Plane",
          plotOutput(outputId ="SO_CE_plane")
        )
      ) # End tabsetPanel
    ) # End mainPanel
  )#close side barlayout
)#close UI fluidpage


server = function(input, output){
  
  
  # This is only used to test the f_wrapper function
  #model_outputs_default = f_wrapper(age_range = "55-64",
  #                                 gender =  "female",  # maximum age of follow up default is 110
  #                                  treatment_names = treatment_names_global[c(1:2)],
  #                                 n_samples = 10)
  
  #when action button pressed ...
  data <- reactive({
    
    treatment_names_local <- unlist(input$SI_treatment_names)
    n_treatments_local <- length(treatment_names_local)
    
    model_outputs = f_wrapper(age_range = input$SI_age_range,
                              gender = input$SI_gender,  # maximum age of follow up default is 110
                              treatment_names = input$SI_treatment_names,
                              n_samples = input$SI_n_samples)
    
    return(list("treatment_names_local" = treatment_names_local,
                "n_treatments_local" = n_treatments_local,
                "model_outputs" = model_outputs))
  })
  
  
  #—— CREATE COST EFFECTIVENESS TABLE ——#
  
  
  #renderTable continuously updates table
  output$SO_icer_table =renderTable({
    
    format_results <- function(x, n_digits = 2) {
      paste0(format(mean(x), nsmall = n_digits, digits = n_digits), " (",
             format(quantile(x, prob = 0.025), nsmall = n_digits, digits = n_digits), " ,",
             format(quantile(x, prob = 0.975), nsmall = n_digits, digits = n_digits), ")")
    }
    
    results_matrix <- matrix(nrow = data()$n_treatments_local, ncol = 4)
    rownames(results_matrix) <- data()$treatment_names_local
    for(treatment in data()$treatment_names_local) {
      results_matrix[treatment, ] <- c(format_results(data()$model_outputs$total_qalys[treatment,]),
                                       format_results(data()$model_outputs$total_costs[treatment,]),
                                       format_results(data()$model_outputs$net_benefit[treatment,]),
                                       format_results(data()$model_outputs$incremental_net_benefit[treatment,]))
    }
    
    
    results_matrix = cbind(results_matrix, round(data()$model_outputs$ICER,2))
    colnames(results_matrix)  = c("QALYs", "Costs", "Net_benefit", "Inc.Net_benefit", "ICER")
    
    
    for (treatment in data()$treatment_names_local){
      if(results_matrix[treatment,5]<0) {
        if (mean(data()$model_outputs$total_qalys[treatment,])-mean(data()$model_outputs$total_qalys[1,])>0){results_matrix[treatment,5] = "Dominant"}
        if (mean(data()$model_outputs$total_qalys[treatment,])-mean(data()$model_outputs$total_qalys[1,])<0){results_matrix[treatment,5] = "Dominated"} 
      }
      
    }
    
    icer_table = data.frame(results_matrix)
    icer_table
    #close data—frame
    #round the data—frame to two digits
    
    #print the results table
    # df_res_table
    
  }, rownames = TRUE) 
  
  
  #table plot end.
  #—— CREATE COST EFFECTIVENESS PLANE ——#
  
  #render plot repeatedly updates.
  output$SO_CE_plane =renderPlot({
    
    #calculate incremental costs and qalys
    knips_bcea <- bcea(e = t(data()$model_outputs$total_qalys), c = t(data()$model_outputs$total_costs), 
                       ref = 1, interventions = data()$treatment_names_local) 
    #create cost effectiveness plane plot
    knips_multi_ce <- multi.ce(knips_bcea)
    ceac.plot(knips_multi_ce, graph = "ggplot",
              line = list(colors = rainbow(12)),
              pos = c(1,1))
  })#renderplot end
  
}

shinyApp(ui , server)

