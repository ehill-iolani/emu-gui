library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(plotly)

# Define UI for application, create sidebar and body
ui <- shinyUI(fluidPage(
  titlePanel("AIN EMU 16S Classifier (beta)"),
  sidebarLayout(
    sidebarPanel(
      helpText("Select fastq or fastq.gz files you want to classify"),
      fileInput("fastqs", "Upload .fastq or fastq.gz files", accept = c(".fastq", ".gz"), multiple = TRUE),
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
  observeEvent(input$run, {
  showModal(modalDialog("Running the EMU pipeline, please be patient...ingesting files", footer = NULL))

    # Make a directory to store input files
    dir.create(file.path("/home/staging"))

    # Set large upload size limit (server side)
    options(shiny.maxRequestSize = 70 * 1024^2)

    # Make a directory for each barcode and
    # copy the fastq files to their respective barcode directories
    files <- input$fastqs$name
    bc_names <- str_extract_all(input$fastqs$name, "barcode[0-9]+") %>% unique()
    bc_names <- unlist(bc_names)
    for (i in seq_along(bc_names)) {
      dir.create(file.path("/home/staging", bc_names[i]))
    }

    for (i in seq_along(bc_names)) {
      bucket <- which(str_detect(files, bc_names[i]))
      file.copy(input$fastqs$datapath[bucket], file.path("/home/staging", bc_names[i]))
    }

    # Check if the files are fastq or fastq.gz; unzip if necessary
    setwd("/home/staging")
    if(any(str_detect(files, ".gz"))) {
      system("gunzip barcode*/*.gz")
      system("for i in */; do cd $i; for j in *; do mv $j ${j}.fastq; done; cd ..; done")
    }

    removeModal()
    showModal(modalDialog("Running the EMU pipeline, please be patient...filtering files", footer = NULL))

    # Create processing directory and set working directory
    # Cat all the fastqs in a directory into one file and
    # write to processing directory
    dir.create(file.path("/home", "processing"))
    system("for i in */; do j=${i%/}; cat $j/*.fastq >> /home/processing/${j}_unfiltered.fastq; done")

    # Filter reads by length and quality using nanofilt
    setwd("/home/processing")
    system(paste("for i in *.fastq; do j=${i%_*}; conda run -n emu NanoFilt -l", input$lowerlength, "--maxlength", input$upperlength, "-q", input$quality, "$i > ${j}_filtered.fastq; done"))
    system("rm *_unfiltered.fastq")

    # Create a directory to store results
    dir.create(file.path("/home", "results"))

    removeModal()
    showModal(modalDialog("Running the EMU pipeline, please be patient...running EMU", footer = NULL))

    # Run EMU
    system(paste("for i in *.fastq; do conda run -n emu emu abundance --db /home/database --keep-counts --output-unclassified --output-dir ../results --threads", input$threads, "$i; done"))

    removeModal()
    showModal(modalDialog("Running the EMU pipeline, please be patient...visualizing the results", footer = NULL))

    # Read in the EMU results
    setwd("/home/results")
    emu_out <- list.files(pattern = "abundance.tsv")

    # Read in and combine the EMU results
    emu_out_all <- data.frame()
    for (i in emu_out) {
      temp <- read.delim(i, header = TRUE, sep = "\t") %>%
              mutate(barcode = str_extract(i, "barcode[0-9]+"))
      emu_out_all <- rbind(emu_out_all, temp)
    }

    # Plot the relative abundance of each taxon by barcode
    output$results <- renderPlotly({
      ggplotly(ggplot(emu_out_all, aes(x = barcode, y = abundance, fill = genus)) +
            geom_bar(stat = "identity", position = "fill") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
            labs(x = "Barcode", y = "Relative abundance (%)", fill = "Genus") +
            ggtitle("Relative abundance of genera per sample"))
    })

    removeModal()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
