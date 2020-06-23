## From https://www.r-bloggers.com/deploying-an-r-shiny-app-with-docker/
## and https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
##

#------------------------------------------------------------
# Prepare R/Shiny with all packages
#------------------------------------------------------------

FROM rocker/shiny:3.6.3 AS base

RUN apt update && apt install -y apt-utils \
    libcurl4-gnutls-dev libnode-dev \
    libssl-dev libxml2-dev libjpeg-dev \
    libgl-dev libglu-dev tk-dev libhdf5-dev \
    libgit2-dev libssh2-1-dev libnetcdf-dev \
    libudunits2-dev libgdal-dev libbz2-dev \
    jags emacs git procps htop python-pip \
    pdftk

## Remove old folder
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny 
RUN rm -fr /omicsplayground

#------------------------------------------------------------
## Pull from GitHub
#------------------------------------------------------------
WORKDIR /omicsplayground
COPY R /omicsplayground/R
COPY shiny /omicsplayground/shiny

#------------------------------------------------------------
## Install R packages that are required
#------------------------------------------------------------
RUN R -e "setwd('R');source('requirements.R')"
RUN R -e "install.packages(c('umap','corrplot','wordcloud'))"
RUN pip install umap-learn

## FROM base AS rbase

#------------------------------------------------------------
## Pull from GitHub
#------------------------------------------------------------
WORKDIR /
RUN rm -fr /omicsplayground
RUN git clone https://github.com/bigomics/omicsplayground.git && chmod -R ugo+rwX /omicsplayground

#------------------------------------------------------------
# Copy local data into Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
COPY data-examples/* /omicsplayground/data/
RUN R -e "setwd('data');source('init.R')"

#------------------------------------------------------------
# Copy further configuration files into the Docker image
#------------------------------------------------------------
COPY docker/shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY docker/shiny-server.sh /usr/bin/shiny-server.sh
RUN chmod ugo+rwX /tmp && chmod o+t /tmp

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
