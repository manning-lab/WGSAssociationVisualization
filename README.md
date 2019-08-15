# WGSAssociationVisualization

## Regional plots from WGS summary statistics
Maintainer: Manning Lab <br>
Version: 0.1

## Description
This application is designed to produce regional plots of a given search range from a bgzipped WGS summary statistics file stored in your Google bucket. The regional plots thus produced can be downloaded as PNG files.

### Data required to run this application
 - WGS summary statistics file in bgzipped (.gz) format residing in a Google bucket to which you have access. <br>
 The file should contain at least the following columns:
   - MarkerName: unique variant identifier (ex: rsID or chr_pos_ref_alt)
   - chr: chromosome
   - pos: variant position
   - pvalue: p-value
   <br>
 - A tabix indexed file in (.gz.tbi) format in the same location.
 
### Output generated by this application
 - Regional plots in PNG format.

### Prerequisites - 
* [Tabix](http://www.htslib.org/doc/tabix.html)
* [Docker](https://www.docker.com/)
* [Git](https://git-scm.com/)

## Execution
### Steps to be followed - 
#### Building and running docker images
Execute in your terminal - 
 - Clone this git repository to your local machine <br> 
 ```
   git clone https://github.com/manning-lab/WGSAssociationVisualization.git
 ```
 - Change your working directory to WGSAssociationVisualization <br>
 ```
   cd WGSAssociationVisualization/
 ```
 - Build the base docker image <br> 
 ```
   docker build -t manninglab/gcloud-htslib-shiny-image . 
 ```
 - Change working directory to manninglab-wgs-visualization <br> 
 ```
   cd manninglab-wgs-visualization/ 
 ```
 - Build the manninglab/visualization_app image <br> 
 ```
   docker build -t manninglab/visualization_app . 
 ```
 - Run the image manninglab/visualization_app in a docker container called visualization_app <br> 
 ```
   docker run --rm --name visualization_app -ti -p 3838:80 manninglab/visualization_app /bin/bash 
 ```
 
 #### Running a demo session
 - You are now in the docker container. To execute the demo: <br>
 ```
   ./demo_script 
 ```
 
 Navigate to your web browser - 
 - Go to http://127.0.0.1:3838/
 
 Give the plot a minute to generate. <br>
 
 - Optionally, change the range you want to search in. 
 - Click Submit to view the plot
 - Clicking the Download button will download the plot as a .png file (for eg. Regional_plot_20:60900000-61100000.png) 
 
 #### Exiting the demo session
Navigate back to your terminal - 
 - Type Ctrl+C on Windows or Command+C on Mac <br>
 
 If you wish to run the app using your data in the Google bucket, follow the next steps. <br>
 
 If not, to exit the docker session:
 ```
   exit
 ```
 
 #### Starting the application from a Docker container
 - You are still in the docker container. To execute the appfiles, start with the Google Cloud configuration setup <br>
 ```
   ./config_script 
 ```
 - Copy the URL displayed on your screen, similar to:
 ```
  spawn ./login_script.sh
  Go to the following link in your browser:

    https://accounts.google.com/o/oauth2/auth?redirect_uri=urn%3Aietf%3Awg%3Aoauth%3A2.0%3Aoob&prompt=select_account&response_type=code&client_id=764086051850-6qr4p6gpi6hn506pt8esuq83di311hur.apps.googleusercontent.com&scope=https%3A%2F%2Fwww.googleapis.com%2Fauth%2Fuserinfo.email+https%3A%2F%2Fwww.googleapis.com%2Fauth%2Fcloud-platform+https%3A%2F%2Fwww.googleapis.com%2Fauth%2Faccounts.reauth&access_type=offline


  Enter verification code:
 ```
 
Navigate to your web browser -
 - Paste the URL onto your browser's address bar
 - Select the account associated with the Google bucket containing your file
 - Click Allow
 - Copy the displayed code
  
Navigate back to your terminal - 
 - Paste the code and press enter - <br>
 ```
   Enter verification code: 4/nQEUyDXiVSWoOHoO3oE1Tj7PPQRMIaLXJOCM-m_vFWOUi3kkfzhP5_S5s
   
 ```
 - Enter (**Note:** you will not be able to see the typed characters)
 ```
 ++
 ```
 - If everything is correct, an access token, similar to the one below, will be displayed. Copy the token, you will need it when you run the app.
 ```
   Access token
   ya29.GltjB01piVLyObuiAYK0gZmShRuTiXtKCS2BBaSIKa5qWKW6U-baZinGarINyB_9K_tW_zKJhBzwoUKNqFruFIQqxYKRrKE5L6bgPXO-kwk8xGUxwjE9eR1iZ4dK
 ```
#### Using the application
Navigate to your web browser - 
 - Go to http://127.0.0.1:3838/
 - Paste your access token in the input box
 - Enter the path to the summary statistics file - 
 ```
   Something like this - 
   gs://fc-91605a4c-df34-4248-b17v-ca123456e59/wgs-summary-stats-file.txt.gz
 ```
 - Enter the columns numbers for the Marker name (a unique variant identifier), chromosome of the variant, position of the variant and p-value of the variant. For eg. for the demo file, the input would be 1, 2, 3 and 9 respectively.
 - Enter the range you want to search in. (For eg. 20:60900000-61100000)
 - Click Submit to view the plot
 - Clicking the Download button will download the plot as a .png file (for eg. Regional_plot_10:112948590-113048589.png) 
 
#### Exiting the docker container
Navigate back to your terminal - 
 - Type Ctrl+C on Windows or Command+C on Mac
 - Exit the docker session:
 ```
   exit
 ```
  


