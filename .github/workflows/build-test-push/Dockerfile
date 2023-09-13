##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM masierom/omicsplayground-base:ub2204_v2
##FROM bigomics/omicsplayground:v2.3.0

#------------------------------------------------------------
# Install any extra (forgotten...) Ubuntu libs
#------------------------------------------------------------
ENV DEBIAN_FRONTEND noninteractive

#RUN sed -i 's/archive.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list

# Set the locale (uncomment UTF8)
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

#------------------------------------------------------------
# Install any extra (forgotten...) R packages since base
#------------------------------------------------------------
WORKDIR /
#RUN pip3 install umap-learn
RUN R -e "devtools::install_github('bigomics/shinyChatR')"
#------------------------------------------------------------
# Install supporting bigomics stuff (as R packages)
#------------------------------------------------------------

WORKDIR /
RUN R -e "remotes::install_github('bigomics/PCSF',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/playbase',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigdash',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigLoaders',dependencies=FALSE)"
RUN R -e "remotes::install_github('rstudio/bslib',dependencies=FALSE)"

#------------------------------------------------------------
# Download fresh code from GitHub (not an R package)
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
RUN echo Downloading $BRANCH branch

RUN rm -fr /omicsplayground  \
    && git clone -b $BRANCH --single-branch --depth 1 \
       https://github.com/mauromiguelm/omicsplayground.git \   
    && chmod -R ugo+rwX /omicsplayground 
## tag docker version

WORKDIR /omicsplayground
#RUN make version
#COPY FEATURES.md CHANGELOG.md VERSION /omicsplayground/


#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; f
RUN make sass

#------------------------------------------------------------
# Expose port and set entry CMD
#------------------------------------------------------------
EXPOSE 3838
##CMD ["R", "-e", "shiny::runApp('components/app/R',port=3838,host='0.0.0.0')"]
CMD exec R -e "shiny::runApp('components/app/R', port=3838, host='0.0.0.0', launch.browser=0)" 2>&1 | tee -a run.log

# add entrypoint to execute when the docker container starts up

RUN echo "Current directory: $(pwd)" && ls .github/workflows/build-test-push 

# add entrypoint to execute when the docker container starts up

RUN cp .github/workflows/build-test-push/entrypoint.sh .
RUN chmod +x /omicsplayground/entrypoint.sh
RUN cat /omicsplayground/entrypoint.sh

WORKDIR /omicsplayground
ENTRYPOINT ["/omicsplayground/entrypoint.sh"]
