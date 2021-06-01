library(shiny)
# Library to create alerts for incorrect input
library(shinyBS)
# Package created by Tim <https://github.com/manning-lab/WGSregionalPlot>
library(WGSregionalPlot)
# Library to parse JSON files to implement gsutil functions
library(rjson)
# Library to implement toggle sidebar functionality
library(shinyjs)
# Sourcing helper functions file
source("./functions.R")

# User interface of the app
ui <- fluidPage(
  # Toggle functionality
  useShinyjs(),
  titlePanel("Manning Lab - WGS Association Visualization"),
  tabPanel(
    "tab",
    div(id = "Sidebar",
        sidebarPanel(width=3,
          # Adding input boxes and buttons
          fluidRow(
            textInput(
              "bucket",
              h4("Google bucket link"),
              value = "",
              placeholder = "gs://fc-91605a4c-df34-4248-b17v-ca123456e59"
            ),
            actionButton("bucketsubmit", "Submit bucket link"),
            selectInput("gspath", h4("Files: "), "", selected =
                          ""),
            textInput("searchrange", h4("Search Range"), value =
                        "20:60900000-61100000")
          ),
          fluidRow(
            h4("Column numbers for:"),
            column(5, textInput(
              "marker",
              "Variant ID",
              value = "",
              placeholder = "1"
            )),
            column(5, textInput(
              "chr", "Chromosome", value = "", placeholder = "2"
            ))
          ),
          fluidRow(column(
            5, textInput("pos", "Position", value = "", placeholder = "3")
          ),
          column(
            5, textInput("pval", "P-value", value = "", placeholder = "9")
          )),
          fluidRow(
            radioButtons(
              "genbuild",
              h4("Genome build"),
              choices =  c("hg19" = "hg19", "hg38" = "hg38"),
              selected = "hg19"
            ),
            textInput(
              "bucket2",
              h4("Google bucket link for the LD file (optional)"),
              value = "",
              placeholder = "gs://fc-91605a4c-df34-4248-b17v-ca123456e59"
            ),
            actionButton("bucketsubmit2", "Submit bucket link"),
            selectInput("ldpath", h4("LD file (optional) "), "", selected =
                          ""),
            textInput(
              "ldref",
              h4("LD reference variant (optional)"),
              value = "",
              placeholder = "20-61000005-A-G"
            ),
            actionButton("submit", "View plot"),
            downloadButton(outputId = "down", label = "Download the plot")
          )
        )),
    mainPanel(
      actionButton("toggleSidebar", "Toggle sidebar"),
      bsAlert("alert"), # Placeholder to display an alert
      plotOutput("plot", width = "100%", height = "500"),
      tableOutput("topvars")
    )
  )
)

# Server-side code for the app
server <- function(input, output, session)
{
  #setwd("/tmp")
  
  # Toggling the sidebar
  observeEvent(input$toggleSidebar, {
    shinyjs::toggle(id = "Sidebar")
  })
  
  # Observing the bucket link, to .gz file, submission
  observeEvent(input$bucketsubmit, {
    # Closing any alerts that maybe displayed
    closeAlert(session, "errorAlert")
    
    # Generating the access token using gcloud to be used by tabix
    write(
      system(
        "gcloud auth application-default print-access-token",
        intern = T
      ),
      "access_token.txt"
    )
    
    # Handling incorrect bucket link provided by user
    # Initializing warning and error to NULL
    E <- NULL
    W <- NULL
    bucket <- isolate(input$bucket)
    
    # Handling case if user enters link without gs://
    if(!startsWith(input$bucket, "gs://"))
    {
      bucket <- paste0("gs://", input$bucket)
    }
    
    bucket_data <- list(
      # Using tryCatch with error and warning handlers
      value = withCallingHandlers(
        tryCatch(
          read.table(pipe(
            # Using gsutil to obtain list of the bucket contents
            paste(
              "gsutil ls", bucket
            )
          )),
          error = function(e){ E <<- "Incorrect bucket link entered"}
        ),
        warning = function(w) {
          W <<- w
        }
      ),
      warning = W,
      error = E
    )
    
    # If gsutil was executed without errors or warnings
    if(is.null(bucket_data$error) && is.null(bucket_data$warning)){
      # Defaulting the LD data bucket to first bucket link
      updateTextInput(session,
                      "bucket2", 
                      value = bucket)
      
      bucket_items <-
        data.frame(name = as.character(bucket_data$value$V1),
                   stringsAsFactors = FALSE)
      
      # Extracting only .gz files who have corresponding .gz.tbi in the same location
      opt_list <- list()
      for (i in bucket_items$name)
      {
        if (grepl(".gz$", i))
        {
          if (paste0(i, ".tbi") %in% bucket_items$name)
          {
            opt_list <-
              append(opt_list, substring(i, sapply(gregexpr("\\/", i), tail, 1) + 1))
          }
        }
      }
      
      # Adding the extracted files to the drop down list (Updation)
      updateSelectInput(session,
                        "gspath",
                        choices = opt_list,
                        selected = "")
      
      # Extracting the .csv, .tsv and .txt files from the bucket (probable LD files)
      exts <- c("csv", "txt", "tsv")
      opt_list <- list()
      for (i in bucket_items$name)
      {
        if (substring(i, sapply(gregexpr("\\.", i), tail, 1) + 1) %in% exts)
        {
          opt_list <-
            append(opt_list, substring(i, sapply(gregexpr("\\/", i), tail, 1) + 1))
        }
      }
      
      # Adding the extracted files to the drop down list (Updation)
      updateSelectInput(session,
                        "ldpath",
                        choices = opt_list,
                        selected = "")
      
    }
    # If gsutil generated an error
    else
    {
      # Creating an alert
      createAlert(
        session,
        "alert",
        "errorAlert",
        title = "Error",
        content = bucket_data$error,
        append = FALSE
      )
      
      # Resetting the drop down lists
      updateSelectInput(session,
                        "gspath",
                        choices = "",
                        selected = "")
      updateSelectInput(session,
                        "ldpath",
                        choices = "",
                        selected = "")
    }
  })
  
  observeEvent(input$bucketsubmit2, {
    # Closing any alerts that maybe displayed
    closeAlert(session, "errorAlert")
    
    # Handling incorrect bucket link provided by user
    # Initializing warning and error to NULL
    E <- NULL
    W <- NULL
    bucket <- input$bucket2
    # Checking if user cleared the second bucket
    if(bucket != "")
    {
      # Handling case if user enters link without gs://
      if(!startsWith(input$bucket2, "gs://"))
      {
        bucket <- paste0("gs://", input$bucket2)
      }
      
      bucket_data <- list(
        # Using tryCatch with error and warning handlers
        value = withCallingHandlers(
          tryCatch(
            read.table(pipe(
              # Using gsutil to obtain list of the bucket contents
              paste(
                "gsutil ls", bucket
              )
            )),
            error = function(e){ E <<- "Incorrect bucket link entered"}
          ),
          warning = function(w) {
            W <<- w
          }
        ),
        warning = W,
        error = E
      )
      
      # If gsutil was executed without errors or warnings
      if(is.null(bucket_data$error) && is.null(bucket_data$warning)){
        bucket_items <-
          data.frame(name = as.character(bucket_data$value$V1),
                     stringsAsFactors = FALSE)
        
        # Extracting the .csv, .tsv and .txt files from the bucket (probable LD files)
        exts <- c("csv", "txt", "tsv")
        opt_list <- list()
        for (i in bucket_items$name)
        {
          if (substring(i, sapply(gregexpr("\\.", i), tail, 1) + 1) %in% exts)
          {
            opt_list <-
              append(opt_list, substring(i, sapply(gregexpr("\\/", i), tail, 1) + 1))
          }
        }
        
        # Adding the extracted files to the drop down list (Updation)
        updateSelectInput(session,
                          "ldpath",
                          choices = opt_list,
                          selected = "")
        
      }
      
      # If gsutil generated an error
      else
      {
        # Creating an alert
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Error",
          content = bucket_data$error,
          append = FALSE
        )
        
        # Resetting the drop down list
        updateSelectInput(session,
                          "ldpath",
                          choices = "",
                          selected = "")
      }
    }
    
    # Resetting the drop down list if user cleared the bucket link box
    else
    {
      updateSelectInput(session,
                        "ldpath",
                        choices = "",
                        selected = "")
      updateTextInput(session,
                      "ldref", 
                      value = "")
    }
  })
  
  # Rendering the plot 
  output$plot <- renderPlot({
    input$submit
    # If no bucket link is provided ie. using demo files as input
    if (isolate(input$bucket) == "")
    {
      # Obtaining the tabix results
      res_list <-
        isolate(
          get_tabix_df(file = "1kg-t2d.all.assoc.aug12.txt.gz ", searchrange = input$searchrange)
       )
      
      # If tabix resulted in a warning
      if (!is.null(res_list$warning))
      {
        # Creating an alert warning
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Warning",
          content = res_list$warning,
          append = FALSE
        )
        # Clearing the preview table
        output$topvars <- renderTable(data.frame())
      }
      
      # If tabix resulted in an error
      else if (!is.null(res_list$error))
      {
        # Creating an alert warning
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Error",
          content = res_list$error,
          append = FALSE
        )
        # Clearing the preview table
        output$topvars <- renderTable(data.frame())
      }
      
      # If tabix was successful
      else
      {
        # Closing any alerts that maybe displayed
        closeAlert(session, "errorAlert")
        
        # Calling function to initialize the variables
        var_list <<- isolate(var_init(res_list$value, input, output))
        
        # Clearing the preview table
        output$topvars <- renderTable(data.frame())
        
        # Sending the initialized variables to makePlot()
        plot <- makePlot(var_list)
        
        # Displaying the preview table
        out_table(var_list)
        
        # Rendering the plot
        plot
      }
    }
    
    # If user has provided a bucket link
    else
    {
      bucket <- isolate(input$bucket)
      
      # Handling case if user enters link without gs://
      if(!startsWith(input$bucket, "gs://"))
      {
        bucket <- paste0("gs://", input$bucket)
      }
      
      # Reading the acces token file contents to be used by tabix
      accesstoken <- readLines("access_token.txt")
      
      # Appending the bucket link to the file name
      gspath <- paste0(bucket, "/", isolate(input$gspath))
      
      # Stitching together the command to be sent for tabix
      # Exporting the access token since tabix is accessing a file in a Google bucket
      command <-
        isolate(
          paste0(
            "export GCS_OAUTH_TOKEN=",
            accesstoken,
            " ; tabix ",
            gspath,
            " ",
            input$searchrange
          )
        )
      
      # Obtaining the tabix results
      res_list <- isolate(get_tabix_df(command = command))
      
      # If tabix resulted in a warning
      if (!is.null(res_list$warning))
      {
        # Creating an alert warning
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Warning",
          content = res_list$warning,
          append = FALSE
        )
        # Clearing the preview table
        output$topvars <- renderTable(data.frame())
      }
      
      # If tabix resulted in an error
      else if (!is.null(res_list$error))
      {
        # Creating an alert error
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Error",
          content = res_list$error,
          append = FALSE
        )
        # Clearing the preview table
        output$topvars <- renderTable(data.frame())
      }
      
      # If tabix was successful
      else
      {
        # Closing any alerts that maybe displayed
        closeAlert(session, "errorAlert")
        
        # Calling function to initialize the variables
        var_list <<- isolate(var_init(res_list$value, input, output))
        
        # Clearing the preview table
        output$topvars <- renderTable(data.frame())
        
        # Sending the initialized variables to makePlot()
        plot <- makePlot(var_list)
        
        # Displaying the preview table
        out_table(var_list)
        
        # Rendering the plot
        plot
      }
    }
  })
  
  # DOwnloading the rendered plot
  output$down <- downloadHandler(
    # Setting the file name
    filename = function()
    {
      paste0("Regional_plot_", input$searchrange, ".png")
    },
    
    # Generating the plot using the latest variable values
    content = function(file)
    {
      # Initiating a png file
      png(file, width = 1000)
      
      # Sending the latest var_list to makePlot
      makePlot(var_list)
      
      # Closing the png file
      dev.off()
    }
  )
}

# Setting the port for the app
options(shiny.port = 3838)

# Running the app using the ui and server functions defined above
shinyApp(ui = ui, server = server)
