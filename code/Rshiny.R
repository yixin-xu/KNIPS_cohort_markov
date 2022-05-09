# install 'shiny' if haven't already.
## install.packages("shiny")  # necessary if you don't already have the function 'shiny' installed.
library(shiny)
## we need the function shiny installed, this loads it from the library.

# source the wrapper function.
source("./wrapper.R")




ui  <- fluidPage (#creates empty page
  
  # title of app
  titlePanel("Knips cohort Markov Model in Shiny"),
  #layout is a sidebar-layout
  sidebarLayout(
    sidebarPanel(#open sidebar panel
      #input type 
      selectInput(inputId ="SI_gender",label ="Gender",
                  choices = c("female", "male"), selected = "female"),
      selectInput(inputId ="SI_age_range",label ="Age range",
                  choices = c("0-55", "55-64","65-74", "75-84", "85-Inf"), selected = "55-64"),
      checkboxGroupInput(inputId ="SI_treatment_names",
                         label = "Treatment name",
                         choices = paste("Implant", c("MoP Cem CR_Fix Mono", "MoP Cem CR_Fix Mod", "MoP Cem CR_Mob Mod",
                                                      "MoP Cem PS_Fix Mod", "MoP Cem PS_Mob Mod",  "MoP Cem Con_Con Mod",
                                                      "MoP Unc CR_Fix Mod", "MoP Unc CR_Mob Mod", "MoP Unc PS_Fix Mod",
                                                      "MoP Hyb CR_Fix Mod", "OX Cem CR_Fix Mod", "OX Cem PS_Fix Mod")),
                         selected = paste("Implant", c("MoP Cem CR_Fix Mono", "MoP Cem CR_Fix Mod", "MoP Cem CR_Mob Mod",
                                                       "MoP Cem PS_Fix Mod", "MoP Cem PS_Mob Mod",  "MoP Cem Con_Con Mod",
                                                       "MoP Unc CR_Fix Mod", "MoP Unc CR_Mob Mod", "MoP Unc PS_Fix Mod",
                                                       "MoP Hyb CR_Fix Mod", "OX Cem CR_Fix Mod", "OX Cem PS_Fix Mod")),
                         inline = TRUE),
      numericInput(inputId ="SI_n_samples",
                   label ="Number of samples",
                   value = 1000,
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
  
  #when action button pressed ...
  observeEvent(input$run_model,
               ignoreNULL = F, {
                 
                 #Run  model function with Shiny inputs
                
                 gender = input$SI_gender
                 age_range = input$SI_age_range
                 
                 initial_age <- substring(input$SI_age_range, 1, 2)
                 finial_age <- substring(input$SI_age_range, 4, 6)
                 
                 if(is.infinite(finial_age)) {
                   # Above age_range[1]
                   starting_age <- 90 
                 }else{
                   starting_age <- ceiling((as.numeric(finial_age)+as.numeric(initial_age))/2)
                 }
                 
                 treatment_names = input$SI_treatment_names
                 n_treatments = length(input$SI_treatment_names)
                 n_samples = input$SI_n_samples
                 
                 input_parameters <- generate_input_parameters(input$SI_n_samples)
                 
                 model_outputs = generate_net_benefit(input_parameters, 
                                                      lambda = 20000)
            
                 
                 
                 #-- CREATE COST EFFECTIVENESS TABLE --#
                 
                 #renderTable continuously updates table
                 output$SO_icer_table =renderTable({
                   
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
                      if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[2,])>0){results_matrix[treatment,6] = "Dominant"}
                      if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[2,])<0){results_matrix[treatment,6] = "Dominated"} 
                    }
                    
                  }
                  
                  icer_table = data.frame(results_matrix)
                     
                    #close data-frame
                   #round the data-frame to two digits
                   
                   #print the results table
                   # df_res_table
                   
                 }) 
                 #table plot end.
                 #-- CREATE COST EFFECTIVENESS PLANE --#
                 
                 #render plot repeatedly updates.
                 output$SO_CE_plane =renderPlot({
                   
                   #calculate incremental costs and qalys
                   knips_bcea <- bcea(e = t(model_outputs$total_qalys), c = t(model_outputs$total_costs), 
                                      ref = 2, interventions = input$SI_treatment_names) 
                   #create cost effectiveness plane plot
                   knips_multi_ce <- multi.ce(knips_bcea)
                   ceac.plot(knips_multi_ce, graph = "ggplot",
                             line = list(colors = rainbow(12)),
                             pos = c(1,1))
                 })#renderplot end
                 
               })#Observe event end
}#Server end

shinyApp(ui , server)

