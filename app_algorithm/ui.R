# ui.R script for OTTER application
# The user-interface and widget input for the Shiny application is defined here
# Sends user-defined input to server.R, calls created output from server.R
# Now using shinydashboard for the user-interface
# ------------------------------------------------------------------------------
## Header: add an OTTER logo here at some point
header <- dashboardHeader(
  title = div("OTTER", img(src = "title.png", width = 80))
)  #dashboardHeader

## Sidebar: Can contain tabs, can also contain inputs if desired
# Determines how many tabs exist in the body of the ui
sidebar <- dashboardSidebar(disable = TRUE)
# sidebar <- dashboardSidebar(
#   sidebarMenu(
#     menuItem("PK Data Entry", tabName = "conc-tab", icon = icon("medkit"))  ###placeholder icon###
#   )  #sidebarMenu
# )  #dashboardSidebar

## UI Chunks
# Chunks of ui code to be referred to in ui body to reduce clutter
uiActionButtons <- function(suffix) {
	div(class = "row-fluid",
		div(class = "MyClass", actionButton(paste0("add", suffix), "Add")),
		tags$head(tags$style(type = "text/css", ".MyClass {display: inline-block}")),
		div(class = "MyClass", actionButton(paste0("rem", suffix), "Remove")),
		tags$head(tags$style(type = "text/css", ".MyClass {display: inline-block}")),
		div(class = "MyClass", actionButton(paste0("save", suffix), "Save")),
		tags$head(tags$style(type = "text/css", ".MyClass {display: inline-block}"))
	)  #div
}

## Body: Main content of each tab, as determined by sidebar
body <- dashboardBody(
  # tabItems(
  # ## First Tab: Times and Concentrations
  #   tabItem(tabName = "conc-tab",
      fluidRow(
        box(status = "primary", title = "Individual or Mean-Pooled Data", width = 5,
          actionButton("console","Debug Console"),
          numericInput("nobs", "Desired Number of Observations", value = 9),
          radioButtons("absorp", "Route of Administration", inline = T,
            choices = list("IV Bolus" = 0, "Oral Dose" = 1)
          ),  #radioButtons
          uiOutput("samptimeui"),
          uiActionButtons("samp")
        ),  #box
        box(status = "primary", title = "Plot", width = 7,
          plotOutput("plot")
        )  #box
      )  #fluidRow
  #   )  #conc-tab
  # )  #tabItems
)  #dashboardBody

#ui end
dashboardPage(header, sidebar, body)