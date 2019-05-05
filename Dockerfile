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
RUN mkdir -p /playground/ext-pkg/
WORKDIR /playground

## Upload some packages/files that are needed to the image
COPY ext-pkg/nclust1_1.9.4.tar.gz \
     ext-pkg/nclust_2.1.1.tar.gz \
     ext-pkg/pathview_1.16.7.tar.gz \
     ext-pkg/FARDEEP_1.0.1.tar.gz \
     ext-pkg/Seurat_v2.3.3.tar.gz \
     ext-pkg/

# Install R packages that are required
COPY require.R /playground/
RUN R -e "source('/playground/require.R')"

# Some packages didn't install correctly...
RUN R -e "install.packages(c('fpc'))"

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
