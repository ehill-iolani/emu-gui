library(shiny)
#library(shinydashboard)
library(ggplot2)
library(dplyr)
library(stringr)
#library(tidyr)
library(plotly)

# Define UI for application, create sidebar and body
ui <- shinyUI(fluidPage(
  titlePanel("AIN EMU 16S Classifier"),
  sidebarLayout(
    sidebarPanel(
      helpText("Select fastq or fastq.gz files you want to classify"),
      fileInput("fastqs", "Upload .fastq or fastq.gz files", accept = c(".fastq", ".fastq.gz"), multiple = TRUE),
      numericInput("lowerlength", "Minimum amplicon length", value = 1300),
      numericInput("upperlength", "Maximum amplicon length", value = 1700),
      numericInput("quality", "Minimum Quality Score", value = 10),
      numericInput("threads", "Threads", value = 2),
      actionButton("run", "Run EMU!")),
    mainPanel(
      plotlyOutput("results"))
  )
))

# Define server logic required to draw a histogram
server <- function(input, output) {
  # When the run button is clicked, execute the workflow
  eventReactive(input$run, {
  showModal(modalDialog("Running the EMU pipeline, please be patient...", footer = NULL))

    # Make a directory to store input files
    dir.create(file.path("/home/staging"))

    # Set large upload size limit (server side)
    options(shiny.maxRequestSize = 70 * 1024^2)

    # Make a directory for each barcode
    files <- input$fastqs$name
    bc_names <- str_extract_all(input$fastqs$name, "barcode[0-9]+") %>% unique()
    bc_names <- unlist(bc_names)
    for (i in seq_along(bc_names)) {
      dir.create(file.path("/home/staging", bc_names[i]))
    }

    # Write fastq files to their respective barcode directories
    for (i in seq_along(bc_names)) {
      bucket <- which(str_detect(files, bc_names[i]))
      file.copy(input$fastqs$datapath[bucket], file.path("/home/staging", bc_names[i]))
    }

    # Set the temporary working directory
    setwd("/home/staging")

    # Check if the files are fastq or fastq.gz; unzip if necessary
    if(any(str_detect(files, ".gz"))) {
      system("gunzip barcode*/*.gz")
    }

    # Create processing directory and set working directory
    dir.create(file.path("/home", "processing"))

    # Cat all the fastqs in a directory into one file and
    # write to processing directory
    system("for i in barcode*; do cat $i/*.fastq >> ../processing/${i}_unfiltered.fastq; done")

    # Set working directory to processing
    setwd("/home/processing")

    # Filter reads by length and quality using nanofilt
    system(paste("for i in barcode*; do NanoFilt -l", input$lowerlength, "-L", input$upperlength, "-q", input$quality, "-t", input$threads, "-o", "${i}_filtered.fastq", "$i.fastq; done"))

    # Remove the unfiltered fastq files
    system("rm *_unfiltered.fastq")

    # Create a directory to store results
    dir.create(file.path("/home", "results"))

    # Run EMU
    system("for i in *.fastq; emu abundance --keep-counts --output-unclassified --output-dir ../results --threads 24 $i; done")

    # Set working directory to results
    setwd("/home/results")

    removeModal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
