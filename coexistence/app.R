# Load necessary libraries
library(shiny)
library(deSolve)
library(ggplot2)

# 1. Define the User Interface (UI)
ui <- fluidPage(
  titlePanel("Lotka-Volterra Competition Model"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Initial Populations"),
      numericInput("N1_0", "Species 1 Initial (N1):", value = 10, min = 0),
      numericInput("N2_0", "Species 2 Initial (N2):", value = 10, min = 0),
      
      h4("Growth Rates"),
      sliderInput("r1", "Species 1 (r1):", min = 0, max = 2, value = 0.5, step = 0.05),
      sliderInput("r2", "Species 2 (r2):", min = 0, max = 2, value = 0.5, step = 0.05),
      
      h4("Carrying Capacities"),
      sliderInput("K1", "Species 1 (K1):", min = 10, max = 500, value = 100, step = 10),
      sliderInput("K2", "Species 2 (K2):", min = 10, max = 500, value = 100, step = 10),
      
      h4("Competition Coefficients"),
      sliderInput("alpha12", "Effect of Sp 2 on Sp 1 (\u03B112):", min = 0, max = 3, value = 0.5, step = 0.1),
      sliderInput("alpha21", "Effect of Sp 1 on Sp 2 (\u03B121):", min = 0, max = 3, value = 0.5, step = 0.1),
      
      h4("Simulation Time"),
      numericInput("t_max", "Max Time:", value = 100, min = 10)
    ),
    
    mainPanel(
      h3("Coexistence Criteria Status"),
      uiOutput("criteria_ui"),
      hr(),
      h3("Population Dynamics Over Time"),
      plotOutput("pop_plot")
    )
  )
)

# 2. Define the Server Logic
server <- function(input, output) {
  
  # Function holding the differential equations for deSolve
  lv_comp <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dN1 <- r1 * N1 * ((K1 - N1 - alpha12 * N2) / K1)
      dN2 <- r2 * N2 * ((K2 - N2 - alpha21 * N1) / K2)
      list(c(dN1, dN2))
    })
  }
  
  # Reactive block that runs the ODE solver whenever a parameter changes
  sim_data <- reactive({
    state <- c(N1 = input$N1_0, N2 = input$N2_0)
    parameters <- c(
      r1 = input$r1, r2 = input$r2,
      K1 = input$K1, K2 = input$K2,
      alpha12 = input$alpha12, alpha21 = input$alpha21
    )
    times <- seq(0, input$t_max, by = 0.5)
    
    # Solve the equations
    out <- ode(y = state, times = times, func = lv_comp, parms = parameters)
    as.data.frame(out)
  })
  
  # Logic to check the stable coexistence criteria mathematically
  output$criteria_ui <- renderUI({
    K1 <- input$K1
    K2 <- input$K2
    a12 <- input$alpha12
    a21 <- input$alpha21
    
    # Prevent division by zero if an alpha is set to 0
    val1 <- ifelse(a21 == 0, Inf, K2 / a21)
    val2 <- ifelse(a12 == 0, Inf, K1 / a12)
    
    cond1 <- K1 < val1
    cond2 <- K2 < val2
    
    # Output different HTML based on whether both criteria are met
    if (cond1 && cond2) {
      status <- h4(style = "color: #2e7d32;", "Stable Coexistence: FULFILLED \u2714")
      msg <- p(sprintf("Both K1 < K2/\u03B121 (%.1f < %.1f) AND K2 < K1/\u03B112 (%.1f < %.1f)", K1, val1, K2, val2))
    } else {
      status <- h4(style = "color: #c62828;", "Stable Coexistence: NOT FULFILLED \u2716")
      msg <- p(sprintf("Conditions are K1 < K2/\u03B121 (%.1f < %.1f) AND K2 < K1/\u03B112 (%.1f < %.1f). At least one is false.", K1, val1, K2, val2))
    }
    
    tagList(status, msg)
  })
  
  # Logic to plot the simulation data
  output$pop_plot <- renderPlot({
    df <- sim_data()
    
    ggplot(df, aes(x = time)) +
      geom_line(aes(y = N1, color = "Species 1"), linewidth = 1.2) +
      geom_line(aes(y = N2, color = "Species 2"), linewidth = 1.2) +
      theme_minimal() +
      labs(x = "Time", y = "Population Size", color = "Legend") +
      scale_color_manual(values = c("Species 1" = "#1976D2", "Species 2" = "#D32F2F")) +
      theme(
        text = element_text(size = 14),
        legend.position = "top"
      )
  })
}

# 3. Run the App
shinyApp(ui = ui, server = server)