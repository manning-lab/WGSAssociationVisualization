library(shiny)
library(shinyBS)
library(WGSregionalPlot)
library(rjson)
library(shinyjs)

get_tabix_df <- function(file = NULL,
                         searchrange = NULL,
                         command = NULL)
{
  W = NULL
  E = NULL
  e.handler <- function(e) {
    cat(file = stderr(), paste(e, collapse = ", "))
    if (!is.null(command))
    {
      comm_list <- unlist(strsplit(command, " "))
      searchrange <- comm_list[length(comm_list)]
    }
    search_list <- unlist(strsplit(searchrange, "[[:punct:]]"))
    
    if (as.numeric(search_list[2]) > as.numeric(search_list[3]))
    {
      E <<-
        "Tabix failed due to reason: Start of searchrange cannot be greater than end of searchrange"
    }
    else if (grep(": no lines available in input", paste(e, collapse = " ")) == 1)
    {
      E <<-
        "Tabix failed due to reason: No values found in given searchrange"
    }
    else
    {
      E <<- paste("Weird error: ", e)
    }
  }
  if (!is.null(command))
  {
    res_list <- list(
      value = withCallingHandlers(
        tryCatch(
          read.table(pipe(command), stringsAsFactors = F),
          error = e.handler
        ),
        warning = function(w) {
          W <<- w
        }
      ),
      warning = W,
      error = E
    )
    write.table(res_list[["value"]], "tabix_out.txt")
    return(res_list)
  }
  else
  {
    res_list <-
      list(
        value = withCallingHandlers(
          tryCatch(
            read.table(pipe(paste("/usr/local/htslib-1.9/bin/tabix", file, searchrange)), stringsAsFactors = F),
            error = e.handler
          ),
          warning = function(w) {
            W <<- w
          }
        ),
        warning = W,
        error = E
      )
    return(res_list)
  }
}

makePlot <- function(temp, input, output, session)
{
  gspath <-
    ifelse(input$gspath == "",
           "1kg-t2d.all.assoc.aug12.txt.gz",
           input$gspath)
  
  ld_data <- NULL
  ldref <- ifelse(input$ldref == "", "20-61000005-A-G", isolate(input$ldref))
  
  if(gspath == "1kg-t2d.all.assoc.aug12.txt.gz")
  {
    ldpath <- "1kg-t2d.chr20_60.9M-61.1M.ld.csv"
    ld_data <- load_ld(file = ldpath,  ld_ref = ldref)
  }
  else if(input$ldpath == "")
  {
    ld_data <- NULL
    ldref <- NULL
  }
  else
  {
    ldpath <- paste0(input$bucket2, "/", isolate(input$ldpath))
    write(system(
      paste(
        "/usr/local/gcloud/google-cloud-sdk/bin/gsutil cat",
        ldpath
      ),
      intern = T
    ), "ldFile.txt")
    ld_data <- load_ld(file = "ldFile.txt",  ld_ref = ldref)
  }
  
  marker <- ifelse(input$marker == "", 1, input$marker)
  chr <- ifelse(input$chr == "", 2, input$chr)
  pos <- ifelse(input$pos == "", 3, input$pos)
  pval <- ifelse(input$pval == "", 9, input$pval)
  
  write.table(temp, "preview1.txt")
  temp[,as.numeric(pval)] <- as.double(as.character(temp[,as.numeric(pval)]))
  #cat(file = stderr(), paste(typeof(temp[2,as.numeric(pval)]), "\n"))
  write.table(temp, "preview2.txt")
  
  search_list <- unlist(strsplit(input$searchrange, "[[:punct:]]"))
  chr_sr <- as.numeric(search_list[1])
  start <- as.numeric(search_list[2])
  end <- as.numeric(search_list[3])
  disp <-
    head(temp[order(temp[, paste0("V", pval)]), c(as.numeric(marker),
                                                  as.numeric(chr),
                                                  as.numeric(pos),
                                                  as.numeric(pval))], 10)
  disp <-
    data.frame(
      "Marker Name" = as.character(disp[, 1]),
      "Chromosome" = as.character(disp[, 2]),
      "Position" = as.character(disp[, 3]),
      "P-value" = disp[, 4]
    )
  colnames(disp) <-
    c("Marker Name", "Chromosome", "Position", "P-value")
  output$topvars <- renderTable(disp, digits = -1)
  
  make_regional_plot(
    chr = chr_sr,
    start = start,
    end = end,
    variant_data = temp,
    variant_chr_column = paste0("V", chr),
    variant_pos_column = paste0("V", pos),
    variant_y_column = paste0("V", pval),
    variant_marker_column = paste0("V", marker),
    genome_build = input$genbuild,
    variant_ld_data = ld_data,
    variant_ld_ref = ldref
  )
}


ui <- fluidPage(
  useShinyjs(),
  titlePanel("Manning Lab - WGS Association Visualization"),
  tabPanel(
    "tab",
    div(id = "Sidebar",
        sidebarPanel(
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
      bsAlert("alert"),
      plotOutput("plot", width = "100%", height = "500"),
      tableOutput("topvars")
    )
  )
)


server <- function(input, output, session)
{
  setwd("/tmp")
  
  observeEvent(input$toggleSidebar, {
    shinyjs::toggle(id = "Sidebar")
  })
  res_list <- get_tabix_df(file = "1kg-t2d.all.assoc.aug12.txt.gz", searchrange = "20:60900000-61100000")
  #assign("res_list", get_tabix_df(file = "1kg-t2d.all.assoc.aug12.txt.gz", searchrange = "20:60900000-61100000"), envir = .GlobalEnv)
  temp <- res_list$value
  output$plot <-
    renderPlot(
      make_regional_plot(
        chr = 20,
        start = 60900000,
        end = 61100000,
        variant_data = temp,
        variant_chr_column = "V2",
        variant_pos_column = "V3",
        variant_y_column = "V9",
        variant_marker_column = "V1",
        genome_build = "hg19",
        variant_ld_data = load_ld(file = "1kg-t2d.chr20_60.9M-61.1M.ld.csv", "20-61000005-A-G"),
        variant_ld_ref = "20-61000005-A-G"
      )
    )
  disp <- head(temp[order(temp[, "V9"]), c(1, 2, 3, 9)], 10)
  disp <-
    data.frame(
      "Marker Name" = as.character(disp[, 1]),
      "Chromosome" = as.character(disp[, 2]),
      "Position" = as.character(disp[, 3]),
      "P-value" = disp[, 4]
    )
  colnames(disp) <-
    c("Marker Name", "Chromosome", "Position", "P-value")
  output$topvars <- renderTable(disp, digits = -1)
  
  observeEvent(input$bucketsubmit, {
    closeAlert(session, "errorAlert")
    Sys.setenv(GOOGLE_APPLICATION_CREDENTIALS = "/tmp/application_default_credentials.json")
    Sys.setenv(BOTO_PATH = "/tmp/.boto")
    js <-
      fromJSON(file = "/tmp/application_default_credentials.json")
    boto_data <-
      paste0(
        "[OAuth2]\nclient_id = ",
        js$client_id,
        "\nclient_secret = ",
        js$client_secret,
        "\n\n[Credentials]\ngs_oauth2_refresh_token = ",
        js$refresh_token
      )
    write(boto_data, file = "/tmp/.boto")
    write(
      system(
        "gcloud auth application-default print-access-token",
        intern = T
      ),
      "access_token.txt"
    )
    E <- NULL
    W <- NULL
    bucket <- isolate(input$bucket)
    if(!startsWith(input$bucket, "gs://"))
    {
      bucket <- paste0("gs://", input$bucket)
    }
    bucket_data <- list(
      value = withCallingHandlers(
        tryCatch(
          read.table(pipe(
            paste(
              "/usr/local/gcloud/google-cloud-sdk/bin/gsutil ls", bucket
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
    if(is.null(bucket_data$error)){
      updateTextInput(session,
                      "bucket2", 
                      value = bucket)
      bucket_items <-
        data.frame(name = as.character(bucket_data$value$V1),
                   stringsAsFactors = FALSE)
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
      updateSelectInput(session,
                        "gspath",
                        choices = opt_list,
                        selected = "")
      
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
      updateSelectInput(session,
                        "ldpath",
                        choices = opt_list,
                        selected = "")
      
    }
    else
    {
      createAlert(
        session,
        "alert",
        "errorAlert",
        title = "Error",
        content = bucket_data$error,
        append = FALSE
      )
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
    closeAlert(session, "errorAlert")
    E <- NULL
    W <- NULL
    bucket <- input$bucket2
    if(bucket != "NULL")
    {
      if(!startsWith(input$bucket2, "gs://"))
      {
        bucket <- paste0("gs://", input$bucket2)
      }
      bucket_data <- list(
        value = withCallingHandlers(
          tryCatch(
            read.table(pipe(
              paste(
                "/usr/local/gcloud/google-cloud-sdk/bin/gsutil ls", bucket
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
      if(is.null(bucket_data$error)){
        bucket_items <-
          data.frame(name = as.character(bucket_data$value$V1),
                     stringsAsFactors = FALSE)
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
        updateSelectInput(session,
                          "ldpath",
                          choices = opt_list,
                          selected = "")
        
      }
      else
      {
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Error",
          content = bucket_data$error,
          append = FALSE
        )
        updateSelectInput(session,
                          "ldpath",
                          choices = "",
                          selected = "")
      }
    }
    else
    {
      updateSelectInput(session,
                        "ldpath",
                        choices = "",
                        selected = "")
    }
  })
  
  output$plot <- renderPlot({
    input$submit
    if (isolate(input$bucket) == "")
    {
      res_list <-
        #isolate(
          get_tabix_df(file = "1kg-t2d.all.assoc.aug12.txt.gz ", searchrange = input$searchrange)
       # )
      if (!is.null(res_list$warning))
      {
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Warning",
          content = res_list$warning,
          append = FALSE
        )
        
        output$topvars <- renderTable(data.frame())
      }
      else if (!is.null(res_list$error))
      {
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Error",
          content = res_list$error,
          append = FALSE
        )
        output$topvars <- renderTable(data.frame())
      }
      else
      {
        closeAlert(session, "errorAlert")
        plot <- makePlot(res_list$value, input, output, session)
        plot
      }
    }
    else
    {
      bucket <- isolate(input$bucket)
      if(!startsWith(input$bucket, "gs://"))
      {
        bucket <- paste0("gs://", input$bucket)
      }
      accesstoken <- readLines("access_token.txt")
      gspath <- paste0(bucket, "/", isolate(input$gspath))
      command <-
        #isolate(
          paste0(
            "export GCS_OAUTH_TOKEN=",
            accesstoken,
            " ; /usr/local/htslib-1.9/bin/tabix ",
            gspath,
            " ",
            input$searchrange
          )
        #)
      cat(file = stderr(), paste("Command:", command))
      res_list <- get_tabix_df(command = command)
      assign("res_list", res_list, envir = .GlobalEnv)
      if (!is.null(res_list$warning))
      {
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Warning",
          content = res_list$warning,
          append = FALSE
        )
        output$topvars <- renderTable(data.frame())
      }
      else if (!is.null(res_list$error))
      {
        createAlert(
          session,
          "alert",
          "errorAlert",
          title = "Error",
          content = res_list$error,
          append = FALSE
        )
        output$topvars <- renderTable(data.frame())
      }
      else
      {
        closeAlert(session, "errorAlert")
        plot <- makePlot(res_list[["value"]], input, output, session)
        plot
      }
    }
  })
  
  output$down <- downloadHandler(
    filename = function()
    {
      paste0("Regional_plot_", input$searchrange, ".png")
    },
    content = function(file)
    {
      png(file, width = 720)
      makePlot(res_list[["value"]], input, output, session)
      dev.off()
    }
  )
}

options(shiny.port = 3838)
shinyApp(ui = ui, server = server)
