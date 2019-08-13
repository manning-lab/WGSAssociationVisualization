library(shiny)

ui <- fluidPage(
		titlePanel("Visualization using R Shiny"),
		sidebarPanel(
			     # Adding input box, buttons and drop-down list
			     fluidRow(
				      textInput("accesstoken", h3("Access token")),
				      textInput("gspath", h3("Enter Path to file"))
				      ),
			     fluidRow(
				      h4("Enter column numbers for -"),
				      column(4, textInput("marker","Markername")),
				      column(4, textInput("chr","Chromosome"))
				      ),
			     fluidRow(
				      column(4, textInput("pos","Position")),
				      column(4, textInput("pval","P-value"))
				      ),
			     fluidRow(
				      textInput("searchrange", h3("Search Range")),
				      actionButton("submit", "Submit"),
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
	observeEvent(input$submit,
	{
		command <- paste0("export GCS_OAUTH_TOKEN=",input$accesstoken," ; /usr/local/htslib-1.9/bin/tabix ", input$gspath," ", input$searchrange)
		#output$selected_opt <- renderText(system(command,intern=TRUE))
		assign("temp", read.table(pipe(command)), envir = .GlobalEnv)
		output$plot <- renderPlot(plot(temp[,as.numeric(input$pos)], -log10(temp[,as.numeric(input$pval)]), xlab="Position", ylab="Negative log of P-value"))
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
