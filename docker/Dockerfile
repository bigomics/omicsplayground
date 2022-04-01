##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground:base
##FROM bigomics/omicsplayground:v2.3.0

#------------------------------------------------------------
# Install any extra (forgotten...) packages since base
#------------------------------------------------------------

# RUN apt install -y xvfb
RUN R -e "remotes::install_github('JohnCoene/firebase@omics')" 
RUN R -e "install.packages(c('shinylogs'))"
RUN R -e "install.packages(c('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz','https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz'))"
RUN R -e "BiocManager::install(c('edgeR','DESeq2'))"
# RUN pip3 install umap-learn

#------------------------------------------------------------
# Download fresh code from GitHub
#------------------------------------------------------------
WORKDIR /
ARG TAG=master
RUN echo Downloading $TAG branch
RUN wget -nv https://github.com/bigomics/omicsplayground/archive/$TAG.zip \
    && unzip $TAG.zip \
    && rm -fr /omicsplayground  \
    && mv omicsplayground-$TAG omicsplayground \
    && chmod -R ugo+rwX /omicsplayground 

#------------------------------------------------------------
# Copy any extra data sets into Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
## RUN R -e "setwd('scripts');source('build-datasets.R')"
## COPY data/geiger2016-arginine.pgx data/
RUN R -e "setwd('data');source('init.R')"

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
COPY docker/shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY docker/shiny-server.sh /usr/bin/shiny-server.sh
RUN chmod ugo+rwX /tmp && chmod o+t /tmp
RUN ln -s /var/log/shiny-server/log /log

EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]
