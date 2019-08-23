#!/bin/bash
set timeout -1

gcloud auth application-default login
read REPLY

cp /root/.config/gcloud/application_default_credentials.json /tmp
chmod a+rx /tmp/application_default_credentials.json

chmod a+x /usr/bin/shiny-server.sh
source /usr/bin/shiny-server.sh
