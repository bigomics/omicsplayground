## From https://www.r-bloggers.com/deploying-an-r-shiny-app-with-docker/
## and https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
##

#------------------------------------------------------------
# Prepare R/Shiny with all packages
#------------------------------------------------------------

FROM rocker/shiny:3.5.2

RUN apt-get update && apt-get install -y apt-utils \
    libcurl4-gnutls-dev libv8-3.14-dev \
    libssl-dev libxml2-dev  libjpeg-dev \
    libgl-dev libglu-dev tk-dev libhdf5-dev \
    libgit2-dev libssh2-1-dev

## Remove old folder, create new
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny
RUN rm -fr /omicsplayground/ && mkdir -p /omicsplayground/
WORKDIR /omicsplayground

## Upload some packages/files that are needed to the image
RUN mkdir -p /omicsplayground/ext/packages/
COPY ext/packages/*.tar.gz ext/packages/

# Install R packages that are required
COPY R /omicsplayground/R
COPY shiny /omicsplayground/shiny
RUN R -e "setwd('R');source('requirements.R')"

#------------------------------------------------------------
# Install all Playground and some data under /omicsplayground
#------------------------------------------------------------
COPY data /omicsplayground/data
COPY lib /omicsplayground/lib
COPY scripts /omicsplayground/scripts

RUN chmod -R ugo+rwX /omicsplayground

#------------------------------------------------------------
# Copy further configuration files into the Docker image
#------------------------------------------------------------
COPY docker/shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY docker/shiny-server.sh /usr/bin/shiny-server.sh

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
