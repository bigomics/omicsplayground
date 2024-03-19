##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground-base:ub2204_v2

#------------------------------------------------------------
# Install any extra (forgotten...) Ubuntu libs
#------------------------------------------------------------
ENV DEBIAN_FRONTEND noninteractive

#RUN sed -i 's/archive.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
# Set the locale (uncomment UTF8)
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

RUN apt update && apt install -y libsodium-dev

#------------------------------------------------------------
# Install any extra (forgotten...) R packages since base
#------------------------------------------------------------
WORKDIR /
COPY docker/Rprofile .Rprofile
#RUN pip3 install umap-learn
RUN R -e "install.packages(c('shinyvalidate','sodium','cookies','ruv'),force=TRUE)"
RUN R -e "remotes::install_version('matrixStats',version='1.1.0',force=TRUE)"
RUN R -e "BiocManager::install(c('GSVA'),force=TRUE)"
RUN R -e "remotes::install_github('renozao/xbioc',dependencies=TRUE)"

#------------------------------------------------------------
# Install supporting bigomics stuff (as R packages)
#------------------------------------------------------------

WORKDIR /
#RUN git clone -b main --single-branch --depth 1 \
#    https://github.com/bigomics/playdata.git \
#    && R CMD INSTALL playdata \
#    && rm -fr playdata

RUN R -e "remotes::install_github('bigomics/shinyChatR',dependencies=FALSE)"  ## rename??
RUN R -e "remotes::install_github('bigomics/playbase',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigdash',dependencies=FALSE)"
RUN R -e "remotes::install_github('the-y-company/bsutils',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/wizardR',dependencies=FALSE)"

#------------------------------------------------------------
# Download fresh code from GitHub (not an R package)
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
RUN echo Downloading $BRANCH branch
RUN rm -fr /omicsplayground  \
    && git clone -b $BRANCH --single-branch --depth 1 \
       https://github.com/bigomics/omicsplayground.git \   
    && chmod -R ugo+rwX /omicsplayground 

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN if test -f .Rprofile; then mv .Rprofile .Rprofile.DISABLED; fi
RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; fi
RUN make sass

#------------------------------------------------------------
# Expose port and set entry CMD
#------------------------------------------------------------
EXPOSE 3838
CMD exec Rscript dev/run_app_headless.R 2>&1 | tee -a run.log
