#!/bin/bash

## wget https://www.shinyproxy.io/downloads/shinyproxy_2.3.1_amd64.deb
## dpkg -i shinyproxy_2.3.1_amd64.deb

## Restart shinyproxy??
##service shinyproxy restart

## Startup orca-server and network with docker-compose
docker-compose -f docker-compose.yml down
nohup docker-compose -f docker-compose.yml up &
sleep 5

## delete instance on port 4000 because shinyproxy takes over on 8080
docker stop playground

