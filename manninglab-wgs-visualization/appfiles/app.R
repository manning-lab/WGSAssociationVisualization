library(shiny)
library(WGSregionalPlot)

makePlot <- function(temp, input, output)
{
	if(startsWith(input$ldpath, "gs:"))
	{
		ld_data <- load_ld(df = input$ldpath, input$ldref)
	}
	else if(input$ldpath != "NULL")
	{
		ld_data <- load_ld(file = input$ldpath, input$ldref)
	}
	else
	{
		ld_data <- NULL
	}

	search_list <- unlist(strsplit(input$searchrange, "[[:punct:]]"))
	chr <- as.numeric(search_list[1])
	start <- as.numeric(search_list[2])
	end <- as.numeric(search_list[3])

	output$topvars <- renderTable(head(temp[order(temp[,paste0("V", input$pval)]),], 10))

	make_regional_plot(chr=chr, start=start, end=end, variant_data=temp, variant_chr_column=paste0("V", input$chr),
			   variant_pos_column=paste0("V", input$pos), variant_y_column=paste0("V", input$pval),
			   variant_marker_column = paste0("V", input$marker),
			   variant_ld_data = ld_data, variant_ld_ref = input$ldref)
}

ui <- fluidPage(
		titlePanel("Visualization using R Shiny"),
		sidebarPanel(
			     # Adding input box, buttons and drop-down list
			     fluidRow(
				      textInput("accesstoken", h3("Access token"), value="NULL"),
				      textInput("gspath", h3("Enter path to file"), value="1kg-t2d.all.assoc.aug12.txt.gz"),
				      textInput("ldpath", h3("Enter path to LD file (optional)"), value="NULL")
				      ),
			     fluidRow(
				      h4("Enter column numbers for -"),
				      column(4, textInput("marker","Markername", value="1")),
				      column(4, textInput("chr","Chromosome", value="2"))
				      ),
			     fluidRow(
				      column(4, textInput("pos","Position", value="3")),
				      column(4, textInput("pval","P-value", value="9"))
				      ),
			     fluidRow(
				      textInput("ldref", h3("LD reference variant"), value="NULL"),
				      textInput("searchrange", h3("Search Range"), value="20:60900000-61100000"),
				      actionButton("submit", "View plot"),
				      downloadButton(outputId="down", label="Download the plot"))
			     ),
		mainPanel(
			  plotOutput("plot"),
			  tableOutput("topvars")
			  )
		)


server <- function(input, output, session) 
{
	setwd("/tmp")
	temp <- read.table(pipe(paste0("/usr/local/htslib-1.9/bin/tabix 1kg-t2d.all.assoc.aug12.txt.gz 20:60900000-61100000"))) 
	output$plot <- renderPlot(make_regional_plot(chr=20, start=60900000, end=61100000, variant_data=temp, variant_chr_column="V2",
						     variant_pos_column="V3", variant_y_column="V9", variant_marker_column = "V1",
			  			     variant_ld_data = NULL, variant_ld_ref = NULL))
	output$topvars <- renderTable(head(temp[order(temp[,"V9"]),], 10))

	observeEvent(input$submit,
	{
		if(input$gspath == "1kg-t2d.all.assoc.aug12.txt.gz")
		{
			temp <- read.table(pipe(paste0("/usr/local/htslib-1.9/bin/tabix ", input$gspath, " ", input$searchrange)))

			output$plot <- renderPlot(makePlot(temp, input, output))
		}
		else
		{
			command <- paste0("export GCS_OAUTH_TOKEN=",input$accesstoken," ; /usr/local/htslib-1.9/bin/tabix ", input$gspath," ", input$searchrange)
			assign("temp", data.frame(read.table(pipe(command))), envir = .GlobalEnv)
			
			output$plot <- renderPlot(makePlot(temp, input, output))
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
					       makePlot(temp, input, output)
					       dev.off()
				       }
				       )	 		
	    }

options(shiny.port = 3838)
shinyApp(ui = ui, server = server)
