###FIRST TIME ONLY:
{
install.packages("DT")
install.packages("shiny")
install.packages("shinydashboard")
install.packages("shinydashboardPlus")
install.packages("ggvis")
install.packages("ggplot2")
install.packages("dplyr")
}
###FIRST TIME ONLY ^^^
{
library(DT)
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(ggvis)
library(ggplot2)
library(dplyr)
}
### Define UI for application that draws a scatter plot with file upload and col select
{shinyUI <- dashboardPage(
  dashboardHeader(title = "Isotopic Peak Allocation"),
  dashboardSidebar(
    width = 350,
    radioButtons(
      "fileType_Input",
      label = h4("Choose File type"),
      choices = list(".csv/txt" = 1),
      selected = 1,
      inline = TRUE
    ),
    fileInput(
      'file1',
      h4('Upload Your Data'),
      accept = c(
        'text/csv',
        'text/comma-separated-values,text/plain',
        '.csv'
      )),
    
    sidebarMenu(
      
      selectInput(inputId = "x", label = "Select x-axis Variable:",
                  choices = NULL),
      
      selectInput(inputId = "y", label = "Select y-axis Variable:",
                  choices = NULL),
      
      
      menuItem("Raw Data", tabName = "RawData",
               icon = icon("table")),
      
      menuItem("Data Summary", tabName = "DataSummary",
               icon = icon("calculator")),
      menuItem("Charts", tabName = "Charts",
               icon = icon("chart-bar"))
      
    )),
  dashboardBody(
    tabItems(
      tabItem( tabName = "RawData",
               fluidRow(
                 box(
                   title = "View Data", 
                   width = NULL,
                   status = "primary", 
                   solidHeader = TRUE,
                   collapsible = TRUE,
                   
                   DT::dataTableOutput('contents'))
                 
               )),
      tabItem(tabName = "DataSummary",
              
              box(verbatimTextOutput("summary"))
              
      ), 
      tabItem(tabName = "Charts",
              box(
                plotOutput('plot',dblclick = "plot_click",width = '700', height = "700", click = "plot1_dblclick",
                           brush = brushOpts(
                             id = "plot1_brush",
                             resetOnNew = TRUE)),
                verbatimTextOutput("info"),
                tableOutput("tibbl"),
                actionButton('save_to_global', "Save Table"))
      )
      
    ))
)
}

### Run the Server with Reactive Functions
{
shinyServer <- function(input, output, session) {
  vals <- reactiveVal(NULL)
  ranges <- reactiveValues(x = NULL, y = NULL)
  # Get the upload file
  myData <- reactive({
    req(input$file1)
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath, header = TRUE)
    data
  })
  
  output$contents <- DT::renderDataTable({
    DT::datatable(myData())       
  })
  
  observe({
    data <- myData()
    updateSelectInput(session, 'x', choices = names(data))
  }) 
  #gets the x variable name, will be used to change the plot legends
  xVarName <- reactive({
    input$x
  }) 
  
  observe({
    data <- myData()
    updateSelectInput(session, 'y', choices = names(data))
  }) 
  #gets the y variable name, will be used to change the plot legends
  yVarName <- reactive({
    input$y
  }) 
  
  output$summary <- renderPrint({
    summary(myData())
  })
  
  #create the table or update the values.
  observeEvent(c(input$plot_click$x, input$plot_click$y), {
    if (is.null(vals())) {
      vals(tibble(x = input$plot_click$x, y = input$plot_click$y))
    } else {
      vals(vals() %>%
           add_row(x = input$plot_click$x, y = input$plot_click$y))
    }
  })
  
  output$tibbl <- renderTable({
    vals()
  })
  
  observeEvent(input$save_to_global, {
    assign(print(input$x), vals(), envir = .GlobalEnv)
  })
  
  output$plot <- renderPlot({
    df <- myData()
    df2 <- vals()
    if (is.null(vals())) {
      ggplot(df, aes_string(x = input$x, y = input$y))+
        geom_line()+
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    } else {
      segment_data = data.frame(
        x = c(vals()$x),
        xend = c(vals()$x), 
        y = numeric(length(vals()$y)),
        yend = c(vals()$y))
      ggplot(df, aes_string(x = input$x, y = input$y))+
        geom_line()+
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
        geom_point(data=df2, aes(x= c(vals()$x), y= c(vals()$y)),color="blue", shape="square")+
        geom_segment(data=segment_data, aes(x=x, y=y, xend=xend, yend=yend), col="blue", lwd=1)
    }
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
    })
}
}

### Run the App
shinyApp(shinyUI, shinyServer)

