##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM masierom/omicsplayground:latest

#------------------------------------------------------------
# Set env variables 
#------------------------------------------------------------

ARG BRANCH=gha_shinytest
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen

#------------------------------------------------------------
# Install de-novo bigomics stuff (as R packages)
#------------------------------------------------------------

WORKDIR /
RUN R -e "remotes::install_github('bigomics/PCSF',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/playbase',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigdash',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigLoaders',dependencies=FALSE)"

#------------------------------------------------------------
# Download fresh code from GitHub (not an R package)
#------------------------------------------------------------
WORKDIR /
RUN echo Downloading $BRANCH branch
RUN wget -nv https://github.com/mauromiguelm/omicsplayground/archive/$BRANCH.zip \
    && rm -fr /omicsplayground*  \
    && unzip $BRANCH.zip \
    && mv omicsplayground-$BRANCH omicsplayground \
    && chmod -R ugo+rwX /omicsplayground \
    && rm $BRANCH.zip   
    
#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; fi
RUN make sass

# add entrypoint to execute when the docker container starts up
RUN cp .github/workflows/build-test-push/entrypoint.sh .
RUN chmod +x /omicsplayground/entrypoint.sh
RUN cat /omicsplayground/entrypoint.sh

WORKDIR /omicsplayground
ENTRYPOINT ["/omicsplayground/entrypoint.sh"]