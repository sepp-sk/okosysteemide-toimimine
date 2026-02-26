library(shiny)
library(ggplot2)
library(reshape2)

# Helper function to shift matrices (simulates periodic 'wrap-around' boundaries)
shift <- function(M, dx, dy) {
  nr <- nrow(M)
  nc <- ncol(M)
  M[(1:nr - 1 - dy) %% nr + 1, (1:nc - 1 - dx) %% nc + 1]
}

# 1. UI Definition
ui <- fluidPage(
  titlePanel("Spatial Rock-Paper-Scissors: Kerr et al. 2002"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Spatial Structure Control"),
      sliderInput("mixing", "Homogenisation Rate (0 = Static Plate, 1 = Well-Mixed Flask):", 
                  min = 0, max = 1, value = 0, step = 0.05),
      helpText("Increasing this value forces the community to interact and disperse globally, destroying local spatial structure."),
      
      hr(),
      h4("Model Parameters (from Paper Box 1)"),
      sliderInput("delta_c", "Death probability of C (\U0394c):", min = 0.1, max = 0.5, value = 0.33, step = 0.01),
      sliderInput("delta_s0", "Base death prob of S (\U0394s,0):", min = 0.1, max = 0.5, value = 0.25, step = 0.01),
      sliderInput("delta_r", "Death probability of R (\U0394r):", min = 0.1, max = 0.5, value = 0.31, step = 0.01),
      sliderInput("tau", "Toxicity of C to S (\U03C4):", min = 0.1, max = 1.0, value = 0.75, step = 0.05),
      
      hr(),
      actionButton("step", "Advance 5 Generations", class = "btn-success", style="width:100%; margin-bottom:10px;"),
      actionButton("reset", "Reset Simulation", class = "btn-danger", style="width:100%;")
    ),
    
    mainPanel(
      fluidRow(
        column(6, plotOutput("gridPlot", height = "400px")),
        column(6, plotOutput("popPlot", height = "400px"))
      ),
      wellPanel(
        h4("Current Dynamics:"),
        textOutput("statusText")
      )
    )
  )
)

# 2. Server Definition
server <- function(input, output, session) {
  
  grid_size <- 75
  
  # Reactive values to store the grid and history
  sim <- reactiveValues(
    grid = matrix(0, grid_size, grid_size),
    history = data.frame(Time = numeric(), C = numeric(), S = numeric(), R = numeric()),
    time = 0
  )
  
  # Initialization function
  init_grid <- function() {
    # 0 = Empty, 1 = C, 2 = S, 3 = R
    new_grid <- matrix(sample(0:3, grid_size * grid_size, replace = TRUE, 
                              prob = c(0.1, 0.3, 0.3, 0.3)), grid_size, grid_size)
    sim$grid <- new_grid
    sim$time <- 0
    sim$history <- data.frame(Time = 0, 
                              C = sum(new_grid == 1), 
                              S = sum(new_grid == 2), 
                              R = sum(new_grid == 3))
  }
  
  # Initialize on startup
  observeEvent(session, { init_grid() }, once = TRUE)
  observeEvent(input$reset, { init_grid() })
  
  # The Cellular Automaton Engine
  observeEvent(input$step, {
    grid <- sim$grid
    
    # Run a few steps per button click to speed up the visual evolution
    for (step in 1:5) {
      is_C <- (grid == 1)
      is_S <- (grid == 2)
      is_R <- (grid == 3)
      
      # 1. Calculate Local Frequencies (Moore Neighborhood)
      count_neighbors <- function(M) {
        shift(M, 1, 0) + shift(M, -1, 0) + shift(M, 0, 1) + shift(M, 0, -1) +
          shift(M, 1, 1) + shift(M, -1, -1) + shift(M, 1, -1) + shift(M, -1, 1)
      }
      
      local_C <- count_neighbors(is_C) / 8
      local_S <- count_neighbors(is_S) / 8
      local_R <- count_neighbors(is_R) / 8
      
      # 2. Calculate Global Frequencies
      total_cells <- grid_size * grid_size
      global_C <- sum(is_C) / total_cells
      global_S <- sum(is_S) / total_cells
      global_R <- sum(is_R) / total_cells
      
      # 3. Blend Local and Global based on Homogenisation Rate
      m <- input$mixing
      f_C <- (1 - m) * local_C + m * global_C
      f_S <- (1 - m) * local_S + m * global_S
      f_R <- (1 - m) * local_R + m * global_R
      
      f_total <- f_C + f_S + f_R
      f_total[f_total == 0] <- 1 # Prevent division by zero
      
      prob_C <- f_C / f_total
      prob_S <- f_S / f_total
      prob_R <- f_R / f_total
      
      # 4. Apply Death Rules
      rand_death <- matrix(runif(total_cells), grid_size, grid_size)
      
      death_C <- is_C & (rand_death < input$delta_c)
      death_R <- is_R & (rand_death < input$delta_r)
      # Toxicity (tau) scales with the density of C in the perceived neighborhood
      death_S <- is_S & (rand_death < (input$delta_s0 + input$tau * f_C))
      
      new_grid <- grid
      new_grid[death_C | death_R | death_S] <- 0
      
      # 5. Apply Reproduction Rules to Empty Cells
      is_empty <- (new_grid == 0)
      rand_repro <- matrix(runif(total_cells), grid_size, grid_size)
      
      col_C <- is_empty & (rand_repro < prob_C)
      col_S <- is_empty & (rand_repro >= prob_C) & (rand_repro < (prob_C + prob_S))
      col_R <- is_empty & (rand_repro >= (prob_C + prob_S)) & (rand_repro < (prob_C + prob_S + prob_R))
      
      new_grid[col_C] <- 1
      new_grid[col_S] <- 2
      new_grid[col_R] <- 3
      
      grid <- new_grid
    }
    
    # Update state
    sim$grid <- grid
    sim$time <- sim$time + 5
    
    # Append to history
    new_row <- data.frame(Time = sim$time, 
                          C = sum(grid == 1), 
                          S = sum(grid == 2), 
                          R = sum(grid == 3))
    sim$history <- rbind(sim$history, new_row)
  })
  
  # Render the Spatial Grid Plot
  output$gridPlot <- renderPlot({
    grid_df <- melt(sim$grid)
    grid_df$value <- factor(grid_df$value, levels = 0:3, labels = c("Empty", "C (Colicin)", "S (Sensitive)", "R (Resistant)"))
    
    ggplot(grid_df, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_manual(values = c("Empty" = "white", "C (Colicin)" = "#d7191c", "S (Sensitive)" = "#2b83ba", "R (Resistant)" = "#abdda4")) +
      theme_void() +
      labs(title = paste("Spatial Lattice (Generation", sim$time, ")"), fill = "Strain") +
      theme(legend.position = "bottom")
  })
  
  # Render the Population Line Plot
  output$popPlot <- renderPlot({
    if (nrow(sim$history) > 0) {
      hist_long <- melt(sim$history, id.vars = "Time", variable.name = "Strain", value.name = "Population")
      ggplot(hist_long, aes(x = Time, y = Population, color = Strain)) +
        geom_line(size = 1.2) +
        scale_color_manual(values = c("C" = "#d7191c", "S" = "#2b83ba", "R" = "#abdda4")) +
        theme_minimal() +
        labs(title = "Population Dynamics", x = "Generations", y = "Total Cells")
    }
  })
  
  # Diagnostic Text
  output$statusText <- renderText({
    if (input$mixing == 0) {
      "Local interactions only (Static Plate). Look for clumps chasing each other across the grid, allowing coexistence."
    } else if (input$mixing == 1) {
      "Fully mixed environment (Flask). The Colicin (C) diffuses everywhere killing Sensitive (S) cells, allowing Resistant (R) cells to dominate."
    } else {
      "Intermediate mixing. Spatial structure is degrading."
    }
  })
}

shinyApp(ui = ui, server = server)