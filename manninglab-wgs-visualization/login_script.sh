#!/bin/bash
set timeout -1

gcloud auth application-default login
read REPLY

cp /root/.config/gcloud/application_default_credentials.json /tmp
chmod a+rx /tmp/application_default_credentials.json

echo "Navigate to http://127.0.0.1:3838/ on your web browser"

chmod a+x /usr/bin/shiny-server.sh
source /usr/bin/shiny-server.sh
