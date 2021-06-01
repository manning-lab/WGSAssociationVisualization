## Helper functions to be used by the app (app.R)
library(shiny)
# Library to create alerts for incorrect input
library(shinyBS)
# Package created by Tim <https://github.com/manning-lab/WGSregionalPlot>
library(WGSregionalPlot)
# Library to parse JSON files to implement gsutil functions
library(rjson)
# Library to implement toggle sidebar functionality
library(shinyjs)

# Function to use tabix command and check for errors/warnings that may arise
get_tabix_df <- function(file = NULL,
                         searchrange = NULL,
                         command = NULL)
{
  # Initializing warning and error to NULL
  W = NULL
  E = NULL
  
  # Error handling function
  e.handler <- function(e) {
    E <<- NULL
    cat(file = stderr(), paste(e, collapse = ", "))
    # Obtaining the searchrange while using command parameter using a cloud file instead of a local file
    if (!is.null(command))
    {
      comm_list <- unlist(strsplit(command, " "))
      searchrange <- comm_list[length(comm_list)]
    }
    search_list <- unlist(strsplit(searchrange, "[[:punct:]]"))
    
    # Error due to start of searchrange being greater than the end
    if (as.numeric(search_list[2]) > as.numeric(search_list[3]))
    {
      E <<-
        "Tabix failed due to reason: Start of searchrange cannot be greater than end of searchrange"
    }
    
    # Error due to no values being found in the given searchrange
    else if (grep(": no lines available in input", paste(e, collapse = " ")) == 1)
    {
      E <<-
        "Tabix failed due to reason: No values found in given searchrange"
    }
    
    # Any other error
    else
    {
      E <<- paste("Exact error:", e)
    }
  }
  
  # using command parameter using a cloud file instead of a local file
  if (!is.null(command))
  {
    res_list <- list(
      # Using tryCatch with error and warning handlers
      value = withCallingHandlers(
        tryCatch(
          read.table(pipe(command), stringsAsFactors = F),
          error = e.handler
        ),
        warning = function(w) {
          W <<- w
        }),
      warning = W,
      error = E
    )
    # Returning a list containing 3 elements - 
    #   1) value (data frame containing tabixed data or error if any)
    #   2) warning (warning line, if any, else NULL)
    #   3) error (error line, if any, else NULL)
    return(res_list)
  }
  
  # using a local file (demo files)
  else
  {
    res_list <-
      list(
        # Using tryCatch with error and warning handlers
        value = withCallingHandlers(
          tryCatch(
            read.table(pipe(
              paste("tabix", file, searchrange)
            ), stringsAsFactors = F),
            error = e.handler
          ),
          warning = function(w) {
            W <<- w
          }),
        warning = W,
        error = E
      )
    # Returning a list containing 3 elements - 
    #   1) value (data frame containing tabixed data or error if any)
    #   2) warning (warning line, if any, else NULL)
    #   3) error (error line, if any, else NULL)
    return(res_list)
  }
}

# Function to:
#   1) Initialize parameters, to be sent to WGSregionalPlot::make_regional_plot
#   2) If no user provided, then default to demo data 
#   3) Call WGSregionalPlot::load_ld to load ld_file if provided
var_init <- function(temp, input, output)
{
  # Google bucket link to the summary statistics file (.gz)
  gspath <-
    ifelse(input$gspath == "",
           "1kg-t2d.all.assoc.aug12.txt.gz",
           input$gspath)
  
  # Handling LD data 
  ld_data <- NULL
  ldref <-
    ifelse(input$ldref == "", "20-61000005-A-G", isolate(input$ldref))
  
  # Setting defaults
  if (gspath == "1kg-t2d.all.assoc.aug12.txt.gz")
  {
    ldpath <- "1kg-t2d.chr20_60.9M-61.1M.ld.csv"
    ld_data <- load_ld(file = ldpath,  ld_ref = ldref)
  }
  
  # If user does not wish to include LD information 
  else if (input$ldpath == "" || input$bucket2 == "")
  {
    ld_data <- NULL
    ldref <- NULL
  }
  
  # Writing the file specified by user to a local file and then loading it using load_ld()
  else
  {
    ldpath <- paste0(input$bucket2, "/", isolate(input$ldpath))
    write(system(
      paste(
        "gsutil cat",
        ldpath
      ),
      intern = T
    ), "ldFile.txt")
    ld_data <- load_ld(file = "ldFile.txt",  ld_ref = ldref)
  }
  
  # Setting local variables based on user input
  marker <- ifelse(input$marker == "", 1, input$marker)
  chr <- ifelse(input$chr == "", 2, input$chr)
  pos <- ifelse(input$pos == "", 3, input$pos)
  pval <- ifelse(input$pval == "", 9, input$pval)
  genbuild <- input$genbuild
  
  # Splitting searchrange into chromosome, start and end values
  search_list <- unlist(strsplit(input$searchrange, "[[:punct:]]"))
  chr_sr <- as.numeric(search_list[1])
  start <- as.numeric(search_list[2])
  end <- as.numeric(search_list[3])
  
  # Creating a list of variables to be used by out_table() and makePlot()
  var_list <- list(temp=temp, gspath=gspath, ld_data=ld_data, ldref=ldref, marker=marker, chr=chr, 
                   pos=pos, pval=pval, genbuild=genbuild, chr_sr=chr_sr, start=start, end=end, output=output)
  return(var_list)
}

# Function to generate a preview table to appear below the plot
out_table <- function(var_list) {
  # Splitting list to global environment variables
  list2env(var_list,envir=.GlobalEnv)
  
  # Taking the top 10 rows with the least p-value
  # Selecting only the marker, chr, pos and pval columns
  disp <-
    head(temp[order(temp[, paste0("V", pval)]), c(as.numeric(marker),
                                                  as.numeric(chr),
                                                  as.numeric(pos),
                                                  as.numeric(pval))], 10)
  
  # Adding LD data column to the preview if provided by the user
  if (!is.null(ld_data))
  {
    # Initializing LD column to NA
    disp$LD <- NA
    
    # Using a temp column called test to make text transformations
    disp$test <- disp[, 1]
    
    # Formating the test column so as to contain only 1 type of punctuation
    disp$test <- gsub("[[:punct:]]", "_", disp$test)
    
    # Checking if either either marker in test or ld_data start with chr and making it uniform
    if (startsWith(disp$test[1], "chr") &
        !startsWith(ld_data$MarkerName[1], "chr"))
    {
      disp$test <- sub("^chr", "", disp$test)
    }
    else if (!startsWith(disp$test[1], "chr") &
             startsWith(ld_data$MarkerName[1], "chr"))
    {
      disp$test <- paste0("chr", disp$test)
    }
    
    # Obtaining matching LD info from the ld_data
    for (i in seq(1, 10)) {
      if (disp$test[i] %in% ld_data$MarkerName)
      {
        disp$LD[i] <- ld_data[ld_data$MarkerName == disp$test[i], "ld"]
      }
      else
        disp$LD[i] <- NA
    }
    
    # Defining the data frame format to be displayed
    disp <-
      data.frame(
        "Marker Name" = as.character(disp[, 1]),
        "Chromosome" = as.character(disp[, 2]),
        "Position" = as.character(disp[, 3]),
        "P-value" = disp[, 4],
        "LD" = as.character(disp$LD)
      )
    
    colnames(disp) <-
      c("Marker Name", "Chromosome", "Position", "P-value", "LD")
  }
  
  # IF LD infor is not provided by the user
  else
  {
    # Defining the data frame format to be displayed without the LD column
    disp <-
      data.frame(
        "Marker Name" = as.character(disp[, 1]),
        "Chromosome" = as.character(disp[, 2]),
        "Position" = as.character(disp[, 3]),
        "P-value" = disp[, 4]
      )
    
    colnames(disp) <-
      c("Marker Name", "Chromosome", "Position", "P-value")
  }
  
  # Displaying the output table
  output$topvars <- renderTable(disp, digits = -1)
}

# Function to call WGSregionalPlot::make_regional_plot to plot the regional plot
makePlot <- function(var_list)
{
  # Splitting list to global environment variables
  list2env(var_list,envir=.GlobalEnv)
  
  # Calling WGSregionalPlot::make_regional_plot and sending the required parameters
  make_regional_plot(
    chr = chr_sr,
    start = start,
    end = end,
    variant_data = temp,
    variant_chr_column = paste0("V", chr),
    variant_pos_column = paste0("V", pos),
    variant_y_column = paste0("V", pval),
    variant_marker_column = paste0("V", marker),
    genome_build = genbuild,
    variant_ld_data = ld_data,
    variant_ld_ref = ldref
  )
}