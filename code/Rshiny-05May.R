
f_wrapper <- function(age_range,gender, treatment_names, n_samples){
  
  treatment_names = as.character(treatment_names)
  gender = as.character(gender)
  age_range = as.character(age_range)
  
  initial_age = substring(age_range, 1, 2)
  finial_age = substring(age_range, 4, 6)
  
  
  if(is.infinite(as.numeric(finial_age))) {starting_age <- 90 }else{
    if (initial_age == "0-") {starting_age <- 53}else{starting_age <- ceiling((as.numeric(finial_age)+as.numeric(initial_age))/2)}
  }
  
  if (initial_age == "0-") {n_cycles <- 50}else{
    n_cycles <- 100- as.numeric(initial_age)}
  
  n_treatments <- length(treatment_names)
  
  n_states <- 8
  state_names <- paste("State", c("Post TKR <3 years", "Post TKR >=3 years < 10 years", "Post TKR >=10 years", 
                                  "Early revision",  "middle revision", "late revision", "second revision", "Death"))
  
  
  
  
  input_parameters <- generate_input_parameters(n_samples, sensitivity = NULL)
  
  # Run the Markov model to get the model outputs
  model_outputs <- generate_net_benefit(input_parameters, 
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
                  choices = list("0-55", "55-64"), 
                  selected = "55-64"),
      checkboxGroupInput(inputId ="SI_treatment_names",
                         label = "Treatment name",
                         choices = list("Implant MoP Cem CR_Fix Mod", "Implant MoP Cem CR_Fix Mono", "Implant MoP Cem CR_Mob Mod"),
                         selected = list("Implant MoP Cem CR_Fix Mod", "Implant MoP Cem CR_Fix Mono", "Implant MoP Cem CR_Mob Mod"),
                         inline = TRUE),
      numericInput(inputId ="SI_n_samples",
                   label ="Number of samples",
                   value = 10,
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
                 
                 model_outputs = f_wrapper(age_range = input$SI_age_range,
                                           gender = input$SI_gender,  # maximum age of follow up default is 110
                                           treatment_names = input$SI_treatment_names,
                                           n_samples = input$SI_n_samples)
                 #model_outputs = f_wrapper(age_range = "55-64",
                 #gender = "female",  # maximum age of follow up default is 110
                 #treatment_names = paste("Implant", c("MoP Cem CR_Fix Mod", "MoP Cem CR_Fix Mono", "MoP Cem CR_Mob Mod")),
                 #n_samples = 10)
                 
                 #—— CREATE COST EFFECTIVENESS TABLE ——#
                 
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
                   
                   
                   results_matrix = cbind(treatment_names,results_matrix, round(model_outputs$ICER,2))
                   colnames(results_matrix)  = c("Treatment_names", "QALYs", "Costs", "Net_benefit", "Inc.Net_benefit", "ICER")
                   
                   
                   for (treatment in treatment_names){
                     if(results_matrix[treatment,6]<0) {
                       if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[1,])>0){results_matrix[treatment,6] = "Dominant"}
                       if (mean(model_outputs$total_qalys[treatment,])-mean(model_outputs$total_qalys[1,])<0){results_matrix[treatment,6] = "Dominated"} 
                     }
                     
                   }
                   
                   icer_table = data.frame(results_matrix)
                   
                   #close data—frame
                   #round the data—frame to two digits
                   
                   #print the results table
                   # df_res_table
                   
                 }) 
                 #table plot end.
                 #—— CREATE COST EFFECTIVENESS PLANE ——#
                 
                 #render plot repeatedly updates.
                 output$SO_CE_plane =renderPlot({
                   
                   #calculate incremental costs and qalys
                   knips_bcea <- bcea(e = t(model_outputs$total_qalys), c = t(model_outputs$total_costs), 
                                      ref = 1, interventions = input$SI_treatment_names) 
                   #create cost effectiveness plane plot
                   knips_multi_ce <- multi.ce(knips_bcea)
                   ceac.plot(knips_multi_ce, graph = "ggplot",
                             line = list(colors = rainbow(12)),
                             pos = c(1,1))
                 })#renderplot end
                 
               })#Observe event end
}#Server end

shinyApp(ui , server)
