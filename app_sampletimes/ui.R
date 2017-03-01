fluidPage(
  h1("Spline application",
    style = "text-align:center"
  ),
  hr(),
  sidebarLayout(
    sidebarPanel(
      selectInput("samp",
        "Sample Times",
        choices = glob.times,
        selected = def.times,
        multiple = TRUE,
        selectize = TRUE
      )
    ),
    mainPanel(
      plotOutput("dataPlot")
    )
  )
)  # fluidPage
