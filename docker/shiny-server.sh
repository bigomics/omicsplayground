#!/bin/sh

# Make sure the directory for individual app logs exists
mkdir -p /var/log/shiny-server
chown shiny.shiny /var/log/shiny-server

## NOTICE we are not using Shiny-server anymore but just R/Shiny
## command. For a single app we do not need Shiny Server. One
## layer less....

##exec shiny-server >> /var/log/shiny-server.log 2>&1
exec R -e "shiny::runApp('/omicsplayground/shiny', port=3838, host='0.0.0.0', launch.browser=FALSE)" >> /var/log/shiny-server/omicsplayground-`date +%Y%m%d-%H%M%S`.log 2>&1
##exec R -e "shiny::runApp('/omicsplayground/shiny', port=3838, host='0.0.0.0', launch.browser=FALSE)" >> /var/log/shiny-server/omicsplayground.log 2>&1
