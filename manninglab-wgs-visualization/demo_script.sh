#!/bin/bash

cp /srv/shiny-server/1kg-t2d.all.assoc.aug12.txt.gz /tmp
cp /srv/shiny-server/1kg-t2d.all.assoc.aug12.txt.gz.tbi /tmp
cp /srv/shiny-server/1kg-t2d.chr20_60.9M-61.1M.ld.csv /tmp

chmod a+x /usr/bin/shiny-server.sh
source /usr/bin/shiny-server.sh
