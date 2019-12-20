## From https://www.r-bloggers.com/deploying-an-r-shiny-app-with-docker/
## and https://www.bjoern-hartmann.de/post/learn-how-to-dockerize-a-shinyapp-in-7-steps/
##

#------------------------------------------------------------
# Prepare R/Shiny with all packages
#------------------------------------------------------------

FROM rocker/shiny:3.6.1

RUN apt-get update && apt-get install -y apt-utils \
    libcurl4-gnutls-dev libv8-3.14-dev \
    libssl-dev libxml2-dev libjpeg-dev \
    libgl-dev libglu-dev tk-dev libhdf5-dev \
    libgit2-dev libssh2-1-dev libnetcdf-dev \
    libudunits2-dev libgdal-dev \
    jags emacs git procps htop python-pip

## Remove old folder
RUN mkdir -p /var/lib/shiny-server/bookmarks/shiny 
RUN rm -fr /omicsplayground

#------------------------------------------------------------
## Pull from GitHub
#------------------------------------------------------------
RUN git clone https://github.com/bigomics/omicsplayground.git && chmod -R ugo+rwX /omicsplayground
WORKDIR /omicsplayground

#------------------------------------------------------------
## Install R packages that are required
#------------------------------------------------------------
RUN pip install umap-learn
RUN R -e "setwd('R');source('requirements.R')"

#------------------------------------------------------------
# Copy local data into Docker image
#------------------------------------------------------------
COPY data/* /omicsplayground/data/
RUN R -e "setwd('scripts');source('update-datasets-info.R')"
RUN R -e "install.packages(c('umap','corrplot','wordcloud'))"

#------------------------------------------------------------
# Copy further configuration files into the Docker image
#------------------------------------------------------------
COPY docker/shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY docker/shiny-server.sh /usr/bin/shiny-server.sh
COPY shiny/resources/rmd_loader.css /usr/local/lib/R/site-library/rmarkdown/rmd/h/rmarkdown/rmd_loader.css
COPY shiny/resources/bigomics-rotating.gif /usr/local/lib/R/site-library/rmarkdown/rmd/h/rmarkdown/rmd_loader.gif 

EXPOSE 3838

CMD ["/usr/bin/shiny-server.sh"]
