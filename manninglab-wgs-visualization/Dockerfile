# Using base image containing installations for gcloud SDK, htslib and shiny server
FROM manninglab/gcloud-htslib-shiny-image 

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY appfiles /srv/shiny-server/

EXPOSE 3838

RUN ln -s /usr/local/gcloud/google-cloud-sdk/bin/gcloud /usr/bin/gcloud

COPY shiny-server.sh /usr/bin/shiny-server.sh

ADD config_script ./
ADD login_script.sh ./

RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]

#CMD ["/usr/bin/shiny-server.sh"]
#ENTRYPOINT ["./config_script"]
