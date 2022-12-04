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
# Install any extra (forgotten...) Ubuntu libs
#------------------------------------------------------------
ENV DEBIAN_FRONTEND noninteractive

# Hirsute (21.04) has been moved..
RUN sed -i 's/archive.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
RUN sed -i 's/security.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list

RUN apt update && apt install -y locales

# Set the locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

#------------------------------------------------------------
# Install any extra (forgotten...) R packages since base
#------------------------------------------------------------

RUN R -e "remotes::install_github('JohnCoene/firebase@omics')" 
RUN R -e "install.packages(c('shinylogs','pryr','shinycssloaders','rpart.plot'))"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz')"
RUN R -e "BiocManager::install(c('edgeR','DESeq2'))"
# RUN pip3 install umap-learn

#------------------------------------------------------------
# Download fresh code from GitHub
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
RUN echo Downloading $BRANCH branch
RUN wget -nv https://github.com/bigomics/omicsplayground/archive/$BRANCH.zip \
    && unzip $BRANCH.zip \
    && rm -fr /omicsplayground  \
    && mv omicsplayground-$BRANCH omicsplayground \
    && chmod -R ugo+rwX /omicsplayground 

#------------------------------------------------------------
# Copy any extra data sets into Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
## RUN R -e "setwd('scripts');source('build-datasets.R')"
## COPY data/geiger2016-arginine.pgx data/
##RUN R -e "setwd('data');source('init.R')"

#------------------------------------------------------------
# Run app
#------------------------------------------------------------
WORKDIR /omicsplayground
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('shiny', port=3838, host='0.0.0.0', launch.browser=0)"]

