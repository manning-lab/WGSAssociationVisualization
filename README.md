# WGSAssociationVisualization

## Prerequisites - 
1) Docker
2) Git

## Steps to be followed - 
Execute in your terminal - 
 1) To clone this repository to your local machine 
 ```git clone ```
 2) docker build -t manninglab/gcloud-htslib-shiny-image .
 3) cd manninglab-wgs-visualization/
 4) docker build -t manninglab/visualization_app .
 5) docker run --name visualization_app -ti -p 3838:80 manninglab/visualization_app /bin/bash
 6) ./install_script
 7) Copy the URL 
 
Navigate to your web browser -
 1) Paste the URL onto your browser's address bar
 2) Select the account associated with the Google Cloud Platform
 3) Click Allow
 4) Copy the displayed code
  
Navigate back to your terminal - 
 1) Paste the code in front of - "Enter verification code:" and press Enter
 2) Type "++" 
 3) If everything is correct, an access token will be displayed. Copy the token, you will need it when you run the app.
  
Navigate to your web browser - 
 1) Type in - http://127.0.0.1:3838/
 2) Paste your access token in the input box
 3) Enter the path to the summary statistics file containing the following columns - 
      MarkerName,chr,pos,ref,alt,minor.allele,maf,mac,n,pvalue,SNPID,BETA,SE,ALTFreq
 4) Enter the range you want to search in. (For eg. 10:112948590-113048589)
 5) Click Submit to view the plot
 6) Click Download to download the plot as a .png file. (For eg. Regional_plot_10:112948590-113048589.png)
 7) You can change the range and click submit to view and download the corresponding plots.
 8) Once you are ready to exit the session, copy the docker links, below the plot, to copy the plot on to your local machine
 
Navigate back to your terminal - 
 1) Type Ctrl+C on Windows or Command+C on Mac
 2) Type exit to exit the docker bash session. Note, your docker container is still running on port 3838.
 3) Paste the copied commands to copy the downloaded plot to your local machine.
  
  


