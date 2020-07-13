#!/bin/sh

# Make sure the directory for individual app logs exists
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

##exec shiny-server >> /var/log/shiny-server.log 2>&1
exec R -e "shiny::runApp('/omicsplayground/shiny', port=3838, \
     host='0.0.0.0', launch.browser=FALSE)" >> \
     /var/log/shiny-server/shiny-playground-`date +%Y%m%d-%H%M%S`.log 2>&1
