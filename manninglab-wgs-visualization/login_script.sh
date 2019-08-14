#!/bin/bash
set timeout -1
cp /srv/shiny-server/1kg-t2d.all.assoc.aug12.txt.gz /tmp
cp /srv/shiny-server/1kg-t2d.all.assoc.aug12.txt.gz.tbi /tmp

gcloud auth application-default login
read REPLY

echo "Access token"
gcloud auth application-default print-access-token

chmod a+x /usr/bin/shiny-server.sh
source /usr/bin/shiny-server.sh
