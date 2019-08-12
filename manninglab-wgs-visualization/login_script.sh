#!/bin/bash

set timeout -1
gcloud auth application-default login
read REPLY

echo "Access token"
gcloud auth application-default print-access-token

chmod a+x /usr/bin/shiny-server.sh
source /usr/bin/shiny-server.sh

ls -al | grep *.png
echo "docker cp [container_id]:/tmp/[filename] filename-CONTAINER_ID.png"
