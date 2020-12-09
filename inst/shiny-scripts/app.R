library(shiny)

#This is a helper function which reads a FASTA file into dataframes
#The code for this function was inspired by the below reference
#reference for this code
#Steipe B., ABC project (.utility 4.07) A Bioinformatics Course: Applied Bioinformatics
#http://steipe.biochemistry.utoronto.ca/abc/index.php/Bioinformatics_Main_Page
readFA <- function(files){
  if (length(files) == 1) { files <- readLines(files) }
  files <- files[! grepl("^$", files)]
  sel <- grep("^>", files)
  myFA <- data.frame(id = files[sel],
                     sequences  = character(length(sel)))

  for (i in seq_along(sel)) {
    first <- sel[i] + 1
    last  <- ifelse(i < length(sel), sel[i + 1] - 1, length(files))
    myFA$sequences[i] <- paste0(files[first:last], collapse = "")
  }
  return(myFA)
}

ui <- fluidPage(

  # App title ----
  titlePanel(tags$h1(tags$b("gscVisualizer:"),"genome difference computing")),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      tags$p("The shiny app allows user to choose a input file, the app
             will calculate the difference among those genes and generatea
             plot shows the distribution for all differences."),

      tags$p("If the input file is genome sequence form, then make selection
             ineach catalog in order to run the shiny app, if the input is
             dot-bracket form, the only two selections you need choose from are
             file type and comparison type"),

      # input
      selectInput("type", "Choose the type of file input:",
                   c("Genome sequences" = "seq",
                     "Dot-bracket form" = "dot")),
      selectInput("seqType","If the input file is sequences, please
                   Choose the sequence type:",
                   c("DNA sequences" = "DNA",
                     "RNA sequences" = "RNA")),
      selectInput("codingType","If the input file is sequences, please
                   Choose the sequence type:",
                   c("Coding sequences" = "coding",
                     "Noncoding sequences" = "noncoding")),
      selectInput("compType","Choose the comparison type:",
                   c("compare to reference" = "cr",
                     "compare every pair" = "cp")),
      fileInput("file1", "Choose FASTA file", accept = ".fa"),

      # actionButton
      actionButton(inputId = "button1", label = "Start!"),

      br(),

    ), # End of side pannel

    mainPanel(

      # Output: Tabet
      tabsetPanel(type = "tabs",
                  tabPanel("Plot and summary",
                           plotOutput("plot"), verbatimTextOutput('textOut')),
                  tabPanel("Data table", verbatimTextOutput('textOut1'))
      )
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {

  #press the start button lead to running this reaction
  processInput <- eventReactive(eventExpr = input$button1, {

    #check file existence
    if (! is.null(input$file1)) {

      inFile <- input$file1

      #If statement for checking use input and calculate the differences,
      #then store the dataframe
      if (input$type == "seq") {
        #input file is in gene sequence form
        if (input$compType == "cr") {
          seqCompareAsFile(input$seqType, input$codingType, inFile$datapath)
        } else if (input$compType == "cp") {
          seqCompareAsFilePair(input$seqType, input$codingType, inFile$datapath)
        }
      } else if (input$type == "dot") {
        #input file is in dot-bracket form
        myFA <- readFA(inFile$datapath)
        if (input$compType == "cr") {
          myFA[["difference"]] <- 0
          for (i in seq_along(myFA$sequences)) {
            dif <- dotComp(myFA$sequences[1], myFA$sequences[i], option = 1)
            myFA$difference[i] <- dif
          }
          myFA
        } else if (input$compType == "cp") {
          myFA[["difference"]] <- 0
          i <- 1
          while (i < length(myFA$sequences)) {
            dif <- dotComp(myFA$sequences[i], myFA$sequences[i + 1], option = 1)
            myFA$difference[i + 1] <- dif
            i <- i + 2
          }
          myFA
        }
      } else {
        ;
      }
    }
  })

  #result plot output
  output$plot <- renderPlot ({
    plotter(processInput())
  })

  #text output
  output$textOut <- renderPrint({
    if (! is.null(processInput()))
      summary(processInput()$difference)
  })

  #Data frame output
  output$textOut1 <- renderPrint({
    if (! is.null(processInput()))
      processInput()
  })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
