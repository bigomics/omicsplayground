## From https://www.r-bloggers.com/deploying-an-r-shiny-app-with-docker/
## and https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
##

#------------------------------------------------------------
# Prepare R/Shiny with all packages
#------------------------------------------------------------

FROM rocker/shiny:3.5.1 

RUN apt-get update && apt-get install -y apt-utils \
    libcurl4-openssl-dev libv8-3.14-dev \
    libssl-dev libxml2-dev  libjpeg-dev \
    libgl-dev libglu-dev tk-dev libhdf5-dev

## ???
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny
RUN mkdir -p /playground/ext/
WORKDIR /playground

## Upload some packages/files that are needed to the image
COPY ext/nclust1_1.9.4.tar.gz \
     ext/nclust_2.1.1.tar.gz \
     ext/pathview_1.16.7.tar.gz \
     ext/FARDEEP_1.0.1.tar.gz \
     ext/Seurat_v2.3.3.tar.gz \
     ext/

# Install R packages that are required
COPY requirements.R /playground/
RUN R -e "source('/playground/requirements.R')"

# Some extra packages so we can use docker cache
COPY requirements2.R /playground/
RUN R -e "source('/playground/requirements2.R')"

#------------------------------------------------------------
# Install all Playground files under /playground
#------------------------------------------------------------

RUN mkdir -p /playground/pgx
COPY pgx /playground/pgx
COPY shiny /playground/shiny
COPY R /playground/R
COPY lib /playground/lib
COPY scripts /playground/scripts

RUN chmod -R ugo+rwX /playground

#------------------------------------------------------------
# Copy further configuration files into the Docker image
#------------------------------------------------------------
COPY docker/shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY docker/shiny-server.sh /usr/bin/shiny-server.sh

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
