#!/bin/bash

cp /srv/shiny-server/1kg-t2d.all.assoc.aug12.txt.gz /tmp
cp /srv/shiny-server/1kg-t2d.all.assoc.aug12.txt.gz.tbi /tmp

chmod a+x /usr/bin/shiny-server.sh
source /usr/bin/shiny-server.sh
