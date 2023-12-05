library(shiny)
library(ggplot2)
library(dplyr)
library(stringr)
library(plotly)
library(shinydashboard)
library(shinyalert)

ui <- dashboardPage(
  dashboardHeader(title = "AIN EMU 16S Classifier (beta)"),
  dashboardSidebar(sidebarMenu(
    menuItem("Emu Classifier & Visualizer", tabName = "emu_class", icon = icon("upload")),
    menuItem("Emu Visualizer", tabName = "emu_viz", icon = icon("upload")),
    menuItem("Relative Abundance", tabName = "relab_plot", icon = icon("chart-simple")),
    menuItem("Bray Curtis PCoA", tabName = "pcoa", icon = icon("chart-simple")),
    menuItem("Rarefaction Curve", tabName = "rare", icon = icon("chart-simple"))
  )),
  dashboardBody(tabItems(
    tabItem(tabName = "emu_class",
            fluidRow(
              box(
                title = "File Upload", status = "primary", solidHeader = TRUE, width = 12,
                fileInput("fastqs", "Upload .fastq or fastq.gz files", accept = c(".fastq", ".gz"), multiple = TRUE),
                  numericInput("lowerlength", "Minimum amplicon length", value = 1300),
                  numericInput("upperlength", "Maximum amplicon length", value = 1700),
                  numericInput("quality", "Minimum Quality Score", value = 10),
                  numericInput("threads", "Threads", value = 2),
                  actionButton("run_class", "Run EMU Classifier + Visualizer!")))),
    tabItem(tabName = "emu_viz",
            fluidRow(
              box(
                title = "File Upload", status = "primary", solidHeader = TRUE, width = 12,
                fileInput("tsv", "Upload EMU .tsv classifications", accept = c(".txt", ".tsv"), multiple = TRUE),
                actionButton("run_viz", "Run EMU Visualizer!")))),
    tabItem(tabName = "relab_plot",
            fluidRow(
              box(
                title = "Relative Abundance",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotlyOutput("relabresults")))),
    tabItem(tabName = "pcoa",
            fluidRow(
              box(
                title = "Bray Curtis PCoA",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotlyOutput("pcoaresults")))),
    tabItem(tabName = "rare",
            fluidRow(
              box(
                title = "Rarefaction Curve",
                status = "primary",
                solidHeader = TRUE,
                width = 12,
                plotlyOutput("rareresults"))))
  ))
)

server <- shinyServer(function(input, output) {
  # Runs the classifier and visualizer
  observeEvent(input$run_class, {
    if (is.null(input$fastqs)) {
      shinyalert(
        title = "Error!",
        text = "Please upload fastq files.",
        type = "error",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
      )
      return(NULL)
    } else {
    showModal(modalDialog("Running the workflow, please be patient...ingesting files", footer = NULL))

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
    showModal(modalDialog("Running the workflow, please be patient...filtering reads using NanoFilt", footer = NULL))

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
    showModal(modalDialog("Running the workflow, please be patient...classifying reads using EMU", footer = NULL))

    # Run EMU
    system(paste("for i in *.fastq; do conda run -n emu emu abundance --db /home/database --keep-counts --output-unclassified --output-dir ../results --threads", input$threads, "$i; done"))

    removeModal()
    showModal(modalDialog("Running the workflow, please be patient...visualizing the using ggplotly", footer = NULL))

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
    output$relabresults <- renderPlotly({
      ggplotly(ggplot(emu_out_all, aes(x = barcode, y = abundance, fill = genus)) +
            geom_bar(stat = "identity", position = "fill") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
            labs(x = "Barcode", y = "Relative abundance (%)", fill = "Genus") +
            ggtitle("Relative abundance of genera per sample"))
    })

    removeModal()
    shinyalert(
      title = "Completed!",
      text = "EMU Classifier + Visualizer has finished running! Results can be viewed under the Relative Abundance, Bray Curtis PCoA, and Rarefaction Curve tabs.",
      type = "success",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
    )
  }})

  # Runs only the visualizer
  observeEvent(input$run_viz, {
    if (is.null(input$tsv)) {
      shinyalert(
        title = "Error!",
        text = "Please upload EMU .tsv files.",
        type = "error",
        showConfirmButton = TRUE,
        confirmButtonText = "OK",
      )
      return(NULL)
    } else {
    showModal(modalDialog("Running the workflow, please be patient...visualizing the using ggplotly", footer = NULL))

    # Read in the EMU classifications and combine into 1 dataframe and parsing out the barcode from filename
    files <- input$tsv
    if (is.list(files)) {
      emu_out_all <- data.frame()
      for (i in seq_along(files$name)) {
        temp <- read.delim(files$datapath[i], header = TRUE, sep = "\t") %>%
                mutate(barcode = str_extract(files$name[i], "barcode[0-9]+"))
        emu_out_all <- rbind(emu_out_all, temp)
      }
    } else {
      emu_out_all <- read.delim(files$datapath, header = TRUE, sep = "\t") %>%
                      mutate(barcode = str_extract(files$name, "barcode[0-9]+"))
    }

    # Plot the relative abundance of each taxon by barcode
    output$relabresults <- renderPlotly({
      ggplotly(ggplot(emu_out_all, aes(x = barcode, y = abundance, fill = genus)) +
            geom_bar(stat = "identity", position = "fill") +
            theme_classic() +
            theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
            labs(x = "Barcode", y = "Relative abundance (%)", fill = "Genus") +
            ggtitle("Relative abundance of genera per sample"))
    })

    removeModal()
    shinyalert(
      title = "Completed!",
      text = "EMU Visualizer has finished running! Results can be viewed under the Relative Abundance, Bray Curtis PCoA, and Rarefaction Curve tabs.",
      type = "success",
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
    )
  }})
})

shinyApp(ui = ui, server = server)
