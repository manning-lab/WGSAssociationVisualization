# Using multi-stage builds feature of Docker (this feature requires Docker 17.05 or higher on the daemon and client)
# Importing and installing Bioconductor packages
FROM bioconductor/release_core2
WORKDIR /tmp
RUN R -e "BiocManager::install(c(\"biomaRt\", \"GenomicRanges\", \"Gviz\", \"EnsDb.Hsapiens.v75\", \"EnsDb.Hsapiens.v86\"), lib=\"/usr/local/lib/R/site-library/\")"
RUN R -e "install.packages(c('data.table', 'devtools', 'draw'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github(\"manning-lab/WGSregionalPlot\")"

FROM rocker/r-ver:3.6.0
# The code below was taken from: https://stackoverflow.com/questions/28372328/how-to-install-the-google-cloud-sdk-in-a-docker-image

RUN apt-get update && apt-get install -qqy \
  curl \
  python-dev \
  python-setuptools \
  wget \
  gcc \
  xtail \
  make \
  libbz2-dev \
  zlib1g-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  gdebi-core \ 
  expect

# Downloading gcloud package
RUN curl https://dl.google.com/dl/cloudsdk/release/google-cloud-sdk.tar.gz > /tmp/google-cloud-sdk.tar.gz

# Installing the package
RUN mkdir -p /usr/local/gcloud \
  && tar -C /usr/local/gcloud -xvf /tmp/google-cloud-sdk.tar.gz \
  && /usr/local/gcloud/google-cloud-sdk/install.sh --quiet 

# Downloading and installing HTSlib
ENV HTSLIB_INSTALL_DIR=/usr/local/htslib-1.9
WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.9.tar.bz2 && \
    cd htslib-1.9/ && \
    ./configure  --enable-plugins --prefix=$HTSLIB_INSTALL_DIR --enable-libcurl --enable-gcs && \
    make && \
    make install && \
    cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/local/ 

# Adding the package path to local
ENV PATH $PATH:/usr/local/gcloud/google-cloud-sdk/bin
ENV PATH $PATH:/usr/local/htslib/htslib-1.9/

# Download and install shiny server
RUN wget --no-verbose https://download3.rstudio.org/ubuntu-14.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb

# Install R packages that are required
# TODO: add further package if you need!
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'shinyBS', 'shinyjs', 'rjson'), repos='http://cran.rstudio.com/')"

# Copying build artifact from the previous stage into this new stage
COPY --from=0 /usr/local/lib/R/site-library/ /usr/local/lib/R/site-library/
