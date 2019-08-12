library(shiny)

ui <- fluidPage(
		titlePanel("Visualization using R Shiny"),
		sidebarPanel(
			     # Adding input box, buttons and drop-down list
			     fluidRow(
				      textInput("accesstoken", h3("Access token")),
				      textInput("gspath", h3("Enter Path to file")),
				      textInput("searchrange", h3("Search Range")),
				      actionButton("submit", "Submit"),
				      actionButton("down", "Download plot")) 
			     ),
		mainPanel(
			  plotOutput("plot"),
			  textOutput("comm"),
			  tableOutput("dwnfiles")
			  )
		)


server <- function(input, output, session) 
{
	setwd("/tmp")
	observeEvent(input$submit,
		     {
			     gspath <- input$gspath
			     accesstoken <- input$accesstoken
			     searchrange <- input$searchrange
			     command <- paste0("export GCS_OAUTH_TOKEN=",accesstoken," ; /usr/local/htslib-1.9/bin/tabix ",gspath," ",searchrange)
			     #output$selected_opt <- renderText(system(command,intern=TRUE))
			     temp <- read.table(pipe(command))
			     # output$selected_opt <- renderTable(temp)
			     
			     output$plot <- renderPlot(plot(temp$V3, -log10(temp$V10), xlab="Position", ylab="Negative log of P-value"))

	 observeEvent(input$down,
		                           {
						   png(paste0("Regional_plot_", searchrange, ".png"))
						   plot(temp$V3, -log10(temp$V10), xlab="Position", ylab="Negative log of P-value")
						   dev.off()
						   output$comm <- renderText("To copy the plot/s to host machine, run the following command/s after exiting the session -")
						   files <- data.frame(file = readLines(pipe("ls -al | grep Regional")))
						   files <- data.frame(Commands = paste0("docker cp visualization_app:/tmp/", substring(files$file, regexpr("Reg", files$file)), " ./", substring(files$file, regexpr("Reg", files$file)), "-CONTAINER_NAME.png"))
						   output$dwnfiles <- renderTable(files)
					   })
	 		})
	 		
	    }

options(shiny.port = 3838)
shinyApp(ui = ui, server = server)
