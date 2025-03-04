if (!interactive()) {
  print("Running from the command line")
  
  # Get the file path from command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    df <- args[1]
  } else {
    stop("Please provide a stats file path as a command-line argument.")
  }
} else {
  print("Running interactively (e.g., in RStudio or R console)")
  # Use a hardcoded path for interactive sessions
  df <- "~/PROBEst/PROBEst/data/param_search/param_stats.csv"
}

# Load required libraries
library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("Interactive gridsearch results for PROBESt v.0.1.2."),
  
  sidebarLayout(
    sidebarPanel(
      textInput("file", "Path to CSV file"),
      uiOutput("group"),
      uiOutput("xvar"),
      uiOutput("yvar"),
      uiOutput("gridbox"),
      uiOutput("gridvar"),
      downloadButton("download_png", "Download as PNG"),
      downloadButton("download_svg", "Download as SVG")
    ),
    
    mainPanel(
      plotlyOutput("plot")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive function to read the uploaded file
  data <- reactive({
    req(input$file)
    read.csv(input$file)

  })
  
  # Dynamically generate UI for selecting x-axis, y-axis, and grid variables
  output$group <- renderUI({
    req(data())
    selectInput("group", "Select colour variable", choices = c(names(data())), selected = "run")
  })
  
  output$gridbox <- renderUI({
    req(data())
    checkboxInput("grid", "Add grid")
  })
  
  output$xvar <- renderUI({
    req(data())
    selectInput("xvar", "Select X-axis variable", choices = names(data()), selected = "iteration")
  })
  
  output$yvar <- renderUI({
    req(data())
    selectInput("yvar", "Select Y-axis variable", choices = names(data()), selected = "mean_hits")
  })
  
  output$gridvar <- renderUI({
    req(data(), input$grid)
    selectInput("gridvar", "Select Grid variable", choices = names(data()), multiple = T)
  })
  
  # Reactive function to create the ggplot
  plot <- reactive({
    req(data(), input$xvar, input$yvar)
    
    p <- ggplot(data(), aes_string(x = input$xvar, y = input$yvar, color = as.factor(input$group), group = input$group)) +
      geom_point() +
      geom_smooth() +
      theme_minimal()
    
    if (input$grid && !is.null(input$gridvar)) {
      p <- p +
        facet_wrap(as.formula(paste("~", paste(input$gridvar, collapse = "+"))), scales = "free")
    }
    
    ggplotly(p)
  })
  
  # Render the plotly plot
  output$plot <- renderPlotly({
    plot()
  })
  
  # Download handlers for PNG and SVG
  output$download_png <- downloadHandler(
    filename = function() {
      "plot.png"
    },
    content = function(file) {
      p <- plot()
      export(p, file)
    }
  )
  
  output$download_svg <- downloadHandler(
    filename = function() {
      "plot.svg"
    },
    content = function(file) {
      p <- plot()
      export(p, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
