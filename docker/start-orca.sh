#!/bin/bash

docker network rm omicsplayground-net
docker network create omicsplayground-net
docker run -d -p 9091:9091 \
       --restart on-failure:10 \
       --network omicsplayground-net \
       --name orca-server \
       quay.io/plotly/orca
