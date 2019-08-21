library(shiny)
library(shinyBS)
library(WGSregionalPlot)

get_tabix_df <- function(file=NULL, searchrange=NULL, command=NULL)
{
  W = NULL
  E = NULL
  if(!is.null(command))
  {
    comm_list <- unlist(strsplit(command, " "))
    res_list <- list(value = withCallingHandlers(tryCatch(read.table(pipe(command)), 
                                                          error=function(e){ E<<-paste("E: ", e) }), 
                                                 warning=function(w){ W<<-w }),warning = W, error = E)
    return(res_list)
  }
  else
  {
    res_list <- list(value = withCallingHandlers(tryCatch(read.table(pipe(paste("/usr/local/htslib-1.9/bin/tabix", file, searchrange))), 
                                                          error=function(e){ E<<-paste("E: ", e) }), 
                                                 warning=function(w){ W<<-w }), warning = W, error = E)
    return(res_list)
  }	
}

makePlot <- function(temp, input, output)
{
  gspath <- ifelse(input$gspath=="", "1kg-t2d.all.assoc.aug12.txt.gz", input$gspath)
  ldpath <- ifelse(input$ldpath=="", "1kg-t2d.chr20_60.9M-61.1M.ld.csv", input$ldpath)
  marker <- ifelse(input$marker=="", 1, input$marker)
  chr <- ifelse(input$chr=="", 2, input$chr)
  pos <- ifelse(input$pos=="", 3, input$pos)
  pval <- ifelse(input$pval=="", 9, input$pval)
  ldref <- ifelse(input$ldref=="", "20-61000005-A-G", input$ldref)
  if(startsWith(ldpath, "gs:"))
  {
    ld_data <- load_ld(df = ldpath, ldref)
  }
  else if(ldpath != "NULL")
  {
    ld_data <- load_ld(file = ldpath, ldref)
  }
  else
  {
    ld_data <- NULL
  }
  
  search_list <- unlist(strsplit(input$searchrange, "[[:punct:]]"))
  chr_sr <- as.numeric(search_list[1])
  start <- as.numeric(search_list[2])
  end <- as.numeric(search_list[3])
  disp <- head(temp[order(temp[,paste0("V", pval)]), c(as.numeric(marker), as.numeric(chr), as.numeric(pos), as.numeric(pval))], 10)
  disp <- data.frame("Marker Name"=as.character(disp[,1]), "Chromosome"=as.character(disp[,2]), "Position"=as.character(disp[,3]), "P-value"=disp[,4])
  colnames(disp) <- c("Marker Name", "Chromosome", "Position", "P-value")
  output$topvars <- renderTable(disp, digits=-1)
  
  make_regional_plot(chr=chr_sr, start=start, end=end, variant_data=temp, variant_chr_column=paste0("V", chr),
                     variant_pos_column=paste0("V", pos), variant_y_column=paste0("V", pval),
                     variant_marker_column = paste0("V", marker),
                     variant_ld_data = ld_data, variant_ld_ref = ldref)
}

ui <- fluidPage(
  titlePanel("Manning Lab - WGS Association Visualization"),
  sidebarPanel(
    # Adding input boxes and buttons
    fluidRow(
      textInput("accesstoken", h3("Access token"), value="", 
                placeholder = "4/nQEUyDXiVSWoOHoO3oE1Tj7PPQRMIaLXJOCM-m_vFWOUi3kkfzhP5_S5s"),
      textInput("gspath", h3("Enter path to file"), value="",
                placeholder = "gs://fc-91605a4c-df34-4248-b17v-ca123456e59/wgs-summary-stats-file.txt.gz"),
      textInput("searchrange", h3("Search Range"), value="20:60900000-61100000")
    ),
    fluidRow(
      h3("Enter column numbers for -"),
      column(4, textInput("marker","Markername", value="", placeholder="1")),
      column(4, textInput("chr","Chromosome", value="", placeholder="2"))
    ),
    fluidRow(
      column(4, textInput("pos","Position", value="", placeholder="3")),
      column(4, textInput("pval","P-value", value="", placeholder="9"))
    ),
    fluidRow(
      textInput("ldpath", h3("Enter path to LD file (optional)"), value="",
                placeholder = "gs://fc-91605a4c-df34-4248-b17v-ca123456e59/ld-data-file.ld.csv"),
      textInput("ldref", h3("LD reference variant (optional)"), value="", placeholder = "20-61000005-A-G"),
      actionButton("submit", "View plot"),
      downloadButton(outputId="down", label="Download the plot"))
  ),
  mainPanel(
    bsAlert("alert"),
    plotOutput("plot"),
    tableOutput("topvars")
  )
)


server <- function(input, output, session) 
{
  setwd("/tmp")
  res_list <- get_tabix_df(file="1kg-t2d.all.assoc.aug12.txt.gz", searchrange="20:60900000-61100000")
  temp <- res_list$value
  output$plot <- renderPlot(make_regional_plot(chr=20, start=60900000, end=61100000, variant_data=temp, 
                                               variant_chr_column = "V2", variant_pos_column = "V3", variant_y_column =
                                                 "V9", variant_marker_column = "V1",
                                               variant_ld_data = load_ld(file = "1kg-t2d.chr20_60.9M-61.1M.ld.csv", "20-61000005-A-G"),
                                               variant_ld_ref = "20-61000005-A-G"))
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
  
  
  output$plot <- renderPlot({
    input$submit
    if(isolate(input$gspath) == "")
    {
      res_list <- isolate(get_tabix_df(file="1kg-t2d.all.assoc.aug12.txt.gz ", searchrange=input$searchrange))
      if(!is.null(res_list$warning))
      {
        createAlert(session, "alert", "errorAlert", title = "Warning", content = res_list$warning, append = FALSE)
        output$topvars <- renderTable(data.frame())
      }
      else if(!is.null(res_list$error))
      {
        createAlert(session, "alert", "errorAlert", title = "Error", content = res_list$error, append = FALSE)
        output$topvars <- renderTable(data.frame())
      }
      else
      {
        closeAlert(session, "errorAlert")
        plot <- isolate(makePlot(res_list$value, input, output))
        plot
      }
    }
    else
    {
      command <- isolate(paste0("export GCS_OAUTH_TOKEN=", isolate(input$accesstoken)," ; /usr/local/htslib-1.9/bin/tabix ", isolate(input$gspath)," ", isolate(input$searchrange)))
      res_list <- isolate(get_tabix_df(command=command))
      if(!is.null(res_list$warning))
      {
        createAlert(session, "alert", "errorAlert", title = "Warning", content = res_list$warning, append = FALSE)
        output$topvars <- renderTable(data.frame())
      }
      else if(!is.null(res_list$error))
      {
        createAlert(session, "alert", "errorAlert", title = "Error", content = res_list$error, append = FALSE)
        output$topvars <- renderTable(data.frame())
      }
      else
      {
        closeAlert(session, "errorAlert")
        plot <- isolate(makePlot(res_list$value, input, output))
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
      png(file)
      makePlot(res_list$value, input, output)
      dev.off()
    }
  )	 		
}

options(shiny.port = 3838)
shinyApp(ui = ui, server = server)
