# Title
# ------------------------------------------------------------------------------
# Non-reactive objects/expressions
# Load package libraries
  library(shiny)
  library(ggplot2)
  library(plyr)
  library(splines)
# Define custom ggplot theme
  theme_bw2 <- theme_set(theme_bw(base_size = 14))

# Specify and simulate using pharmacokinetic model
  CL <- 10
  V <- 50
  KA <- 0.5
  ERR <- 0.3
  dose <- 50
  times <- seq(from = 0,to = 24,by = 0.25)
  conc <- dose*KA/(V*(KA - (CL/V))) * (exp(-CL/V*times) - exp(-KA*times)) *
    (1 + rnorm(n = length(times), mean = 0, sd = ERR))

# ------------------------------------------------------------------------------
# Reactive objects/expressions
  shinyServer(function(input, output) {

  # Dynamic UI for entering knots
    output$Rknots <- renderUI({
      llply(seq_len(input$nknots), function(i) {
        div(
          numericInput(
            inputId = paste0("knot", i),
            label = switch((i == 1) + 1, NULL, "Breakpoints Locations"),
            value = 2 + i*3,
            min = 1,
            max = 23,
            step = 1,
          )  # numericInput
        )  # div
      })  # llply
    })  # output$Rknots

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  # Create data to be used for plotting
    Rpoints <- reactive({
    # Sample from simulated concentrations
      sample.times <- c(0, 0.25, 0.5, 1, 2, 4, 8, 12, 24)  # to be made input
      sample.conc <- conc[times %in% sample.times]
      data.frame(
        time = sample.times,
        conc = sample.conc
      )  # data.frame
    })  # Rpoints()

    Rlines <- reactive({
    # Assign values from outside function
      K <- unlist(
        llply(seq_len(input$nknots), function(i) {
          get("input")[[paste0("knot", i)]]
        })
      )
      deg <- input$degree
      data <- Rpoints()
    # Create spline model
      mod <- lm(conc ~ bs(time, knots = c(K), degree = deg), data = Rpoints())
      test <- seq(0, 24, by = 0.1)  # test model at these times
      test.df <- data.frame(time = test)  # data.frame(times)
      data.frame(
        time = test,
        conc = predict(mod, newdata = test.df)
      )  # data.frame
    })  # Rlines()

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  # Create plot
    output$dataPlot <- renderPlot({
    # Create the basic ggplot2 object for output
      plotobj <- ggplot()
      plotobj <- plotobj + geom_point(aes(x = time, y = conc), data = Rpoints())
      plotobj <- plotobj + geom_line(aes(x = time, y = conc), data = Rlines(),
        size = 1, colour = "red")
    # Print final object
      plotobj
    })  # renderPlot

  })  # shinyServer
