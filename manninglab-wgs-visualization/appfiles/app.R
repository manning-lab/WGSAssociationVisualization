library(shiny)
library(WGSregionalPlot)

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
			  			     variant_ld_data = load_ld(file = "1kg-t2d.chr20_60.9M-61.1M.ld.csv", "20-61000005-A-G"), 
						     variant_ld_ref = "20-61000005-A-G"))
	disp <- head(temp[order(temp[,"V9"]), c(1,2,3,9)], 10)
	disp <- data.frame("Marker Name"=as.character(disp[,1]), "Chromosome"=as.character(disp[,2]), "Position"=as.character(disp[,3]), "P-value"=disp[,4])
	colnames(disp) <- c("Marker Name", "Chromosome", "Position", "P-value")
	output$topvars <- renderTable(disp, digits=-1)

	observeEvent(input$submit,
	{
		if(input$gspath == "")
		{
			temp <- read.table(pipe(paste0("/usr/local/htslib-1.9/bin/tabix 1kg-t2d.all.assoc.aug12.txt.gz ", input$searchrange)))

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
