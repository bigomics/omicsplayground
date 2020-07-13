#!/bin/bash

echo Stopping testing container
docker stop testing && docker rm testing
sleep 3

echo Creating testing container
##docker container create --name testing bigomics/omicsplayground:base
##docker container create --name testing bigomics/omicsplayground:v2.3.0
docker container create --name testing bigomics/omicsplayground:v2.4
docker start testing
sleep 1
##docker exec testing rm -f /omicsplayground/shiny/modules/UsersMapModule.R

echo Deleting old files in container
docker exec testing bash -c 'rm -fr /omicsplayground/*'
docker exec testing bash -c 'rm -f /var/log/shiny-server.log'
docker exec testing bash -c 'rm -f /var/log/shiny-server/*.log'
docker exec testing bash -c 'mkdir -p /root/ShinyApps/log'
docker exec testing bash -c 'mkdir -p /home/shiny/ShinyApps/log'
docker exec testing bash -c 'chmod ugo+rwX /var/log/shiny-server'

echo Copying current git files to container
tar -c `git ls-files` | docker container cp - testing:/omicsplayground
docker exec testing cp docker/shiny-server.conf /etc/shiny-server/shiny-server.conf
#docker exec testing cp docker/shiny-playground.sh /usr/bin/shiny-playground.sh
docker container cp docker/shiny-playground.sh testing:/usr/bin/shiny-playground.sh

PGX="geiger2016-arginine.pgx GSE72056-scmelanoma.pgx"
for p in $PGX; do \
    echo Copying data/$p;    
    docker container cp data/$p testing:/omicsplayground/data; \
done
docker exec testing R -e "setwd('data');source('init.R')"
    
echo Installing extra packages
#docker exec testing cp ext/bin/orca /usr/local/bin
## docker exec testing R -e "install.packages('uwot')"

echo Committing to bigomics/omicsplayground:testing
docker commit testing bigomics/omicsplayground:testing
docker stop testing && docker rm testing

#echo Starting new container at port 4444
#docker run --rm -d -p 4444:3838 --name testing --network docker_play bigomics/omicsplayground:testing


