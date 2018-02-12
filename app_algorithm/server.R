# server.R script for OTTER application
# Reactive objects (i.e., those dependent on widget input) are written here
# ------------------------------------------------------------------------------
server <- function(input, output, session) {

## Setup reactiveValues
# n - counter to determine number of boxes for reactive ui "time, conc" input
# io - on/off switch for processes to be used with data
# df - reactive data.frame that saves values form the "time, conc" input
  rv <- reactiveValues(
    n = 1, io = F,
    df = data.frame(
      time = double(1),
      conc = double(1)
    )  #data.frame
  )  #reactiveValues

## Observe Events
# When pressing add: check values in df are numbers
# If values are numbers then remove default renderUI values and increase values
# Also update dataframe for input - reference values for time, conc input
  observeEvent(input$addsamp, {
    rv$io <- FALSE
    rv$n <- rv$n + 1
    rv$df <- ldply(seq_len(rv$n), function(i) {  #save current values to edit.df
      data.frame(
        time = as.numeric(get("input")[[paste0("time", i)]]),
        conc = as.numeric(get("input")[[paste0("conc", i)]])
      )  #data.frame
    })  #ldply
  })  #observeEvent

# When pressing remove: check to see if > 1 input box, if so decrease values
# Also update dataframe for input - reference values for time,dose,conc input
  observeEvent(input$remsamp, {
    rv$io <- FALSE
    if (rv$n>1) {  #
      rv$n <- rv$n - 1
    }  #if
    rv$df <- ldply(seq_len(rv$n), function(i) {  #save current values to edit.df
      data.frame(
        time = as.numeric(get("input")[[paste0("time", i)]]),
        conc = as.numeric(get("input")[[paste0("conc", i)]])
      )  #data.frame
    })  #ldply
  })  #observeEvent

# When pressing save: update dataframe for input
# Also print to console the state of the data (for debug purposes)
  observeEvent(input$savesamp, {
    rv$df <- ldply(seq_len(rv$n), function(i){  #save current values to edit.df
      data.frame(
        time = as.numeric(get("input")[[paste0("time", i)]]),
        conc = as.numeric(get("input")[[paste0("conc", i)]])
      )  #data.frame
    })  #ldply
    rv$df$time[is.na(rv$df$time)] <- 0
    rv$df$conc[is.na(rv$df$conc)] <- 0
    print(rv$df)
    rv$io <- TRUE
  })  #observeEvent

## UI code for renderUI
# Shows "time, conc" ui as dictated by "n"
# Set up to allow saving of values as they are put in to prevent deletion of 
# progress when Add, Remove is pressed
  inputBoxes <- reactive({
    n <- rv$n
    df <- rv$df
    if(n > 0) {
      llply(seq_len(n), function(i) {
        div(class = "row-fluid",
        # Time
          div(class = "MyClass",  #TIME
            textInput(
              paste0("time", i),
              ifelse(i == 1, "Time", NA),
              df[i, 1]
            )  #textInput
          ),  #input$time-i
          tags$head(tags$style(type = "text/css", 
            ".MyClass {display: inline-block}"
          )),  #make inline
          tags$head(tags$style(type = "text/css", 
            paste0("#time", i, " {max-width: 110px}")
          )),  #change width
        # Conc
          div(class = "MyClass",
            textInput(
              paste0("conc", i),
              ifelse(i == 1, "Concentrations", NA),
              df[i, 2]
            )  #textInput
          ),  #input$conc-i
          tags$head(tags$style(type = "text/css", 
            ".MyClass {display: inline-block}"
          )),  #make inline
          tags$head(tags$style(type = "text/css", 
            paste0("#conc", i, " {max-width: 110px}")
          ))  #change width
        )  #div
      })  #llply
    }  #if
  })  #textboxes

# Finally render the dataframe input UI so it can be used in ui body
  output$samptimeui <- renderUI({inputBoxes()})
  
# Reactive sumexp
  # Rsumexp <- reactive({
  #   
  # })  #Rsumexp
  
# Reactive plot
# Shows data points when rv$io == TRUE
# Shows sumexp curve once available
  # Rplot <- reactive({
  #   
  # })  #Rplot

# Close the R session when Chrome closes
  session$onSessionEnded(function() {
    stopApp()
  })  #endsession

}  #server