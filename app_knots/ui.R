fluidPage(
  h1("Spline application",
    style = "text-align:center"
  ),
  hr(),
  sidebarLayout(
    sidebarPanel(
      sliderInput("degree", "Degree of Piecewise Polynomial",
        value = 2,
        min = 1,
        max = 4
      ),
      numericInput("nknots", "Breakpoints in Spline",
        value = 1,
        min = 1,
        max = 5,
        step = 1
      ),
      uiOutput("Rknots")
    ),
    mainPanel(
      plotOutput("dataPlot")
    )
  )
)  # fluidPage
