library(shiny)

ui <- fluidPage(
		titlePanel("Visualization using R Shiny"),
		sidebarPanel(
			     # Adding input box, buttons and drop-down list
			     fluidRow(
				      textInput("accesstoken", h3("Access token")),
				      textInput("gspath", h3("Enter Path to file"), value="1kg-t2d.all.assoc.aug12.txt.gz")
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
				      textInput("searchrange", h3("Search Range"), value="20:60900000-61100000"),
				      actionButton("submit", "View plot"),
				      downloadButton(outputId="down", label="Download the plot"))
			     ),
		mainPanel(
			  plotOutput("plot")
			  )
		)


server <- function(input, output, session) 
{
	setwd("/tmp")
#	temp <- table(NA)
#	path <- reactive({input$gspath})
#	range <- reactive({input$searchrange})
#	temp <- read.table(pipe(paste0("/usr/local/htslib-1.9/bin/tabix ", path(), " ", range()))) 
#	output$plot <- renderPlot(plot(temp[,as.numeric(input$pos)], -log10(temp[,as.numeric(input$pval)]), xlab="Position", ylab="Negative log of P-value"))
	observeEvent(input$submit,
	{
		if(input$gspath == "1kg-t2d.all.assoc.aug12.txt.gz")
		{
			temp <- read.table(pipe(paste0("/usr/local/htslib-1.9/bin/tabix ", input$gspath, " ", input$searchrange)))
			output$plot <- renderPlot(plot(temp[,as.numeric(input$pos)], -log10(temp[,as.numeric(input$pval)]), xlab="Position", ylab="Negative log of P-value"))
		}
		else
		{
			command <- paste0("export GCS_OAUTH_TOKEN=",input$accesstoken," ; /usr/local/htslib-1.9/bin/tabix ", input$gspath," ", input$searchrange)
			#output$selected_opt <- renderText(system(command,intern=TRUE))
			assign("temp", read.table(pipe(command)), envir = .GlobalEnv)
			output$plot <- renderPlot(plot(temp[,as.numeric(input$pos)], -log10(temp[,as.numeric(input$pval)]), xlab="Position", ylab="Negative log of P-value"))		
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
					       plot(temp[,as.numeric(input$pos)], -log10(temp[,as.numeric(input$pval)]), xlab="Position", ylab="Negative log of P-value")
					       dev.off()
				       }
				       )	 		
	    }

options(shiny.port = 3838)
shinyApp(ui = ui, server = server)
