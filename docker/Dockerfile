##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground-base:ub2204
##FROM bigomics/omicsplayground:v2.3.0

#------------------------------------------------------------
# Install any extra (forgotten...) Ubuntu libs
#------------------------------------------------------------
ENV DEBIAN_FRONTEND noninteractive

#RUN sed -i 's/archive.ubuntu.com/old-releases.ubuntu.com/g' /etc/apt/sources.list
RUN apt update && apt install -y locales \
    librsvg2-dev

# Set the locale (uncomment UTF8)
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

#------------------------------------------------------------
# Install any extra (forgotten...) R packages since base
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN R -e "devtools::install_github('m-jahn/fluctuator')"
RUN R -e "BiocManager::install(c('SBGNview','org.Hs.eg.db','DeMixT'))"
RUN R -e "install.packages(c('svgPanZoom'))"
RUN R -e "install.packages('extdata/packages/pathview_1.16.7.tar.gz')"
# RUN pip3 install umap-learn

#------------------------------------------------------------
# Install bigomics stuff (as R packages)
#------------------------------------------------------------

WORKDIR /omicsplayground
RUN R -e "remotes::install_github('bigomics/bigdash',force=TRUE)"
RUN R -e "devtools::install_github('bigomics/bigLoaders',force=TRUE)"
RUN R -e "devtools::install_github('bigomics/playdata',force=TRUE)"
RUN R -e "devtools::install_github('bigomics/playbase',force=TRUE)"

#------------------------------------------------------------
# Install fresh playbase
#------------------------------------------------------------
# ##RUN R -e "remotes::install_github('bigomics/playbase',force=TRUE)"
# WORKDIR /
# RUN echo Downloading Playbase
# RUN wget -nv https://github.com/bigomics/playbase/archive/main.zip \
#     && rm -fr /playbase-main  \
#     && unzip main.zip \
#     && R CMD INSTALL playbase-main

#------------------------------------------------------------
# Download fresh code from GitHub (not an R package)
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
RUN echo Downloading $BRANCH branch
RUN wget -nv https://github.com/bigomics/omicsplayground/archive/$BRANCH.zip \
    && rm -fr /omicsplayground*  \
    && unzip $BRANCH.zip \
    && mv omicsplayground-$BRANCH omicsplayground \
    && chmod -R ugo+rwX /omicsplayground \
    && rm $BRANCH.zip   

## tag docker version
WORKDIR /omicsplayground
RUN make tag.version

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN if test -f .Rprofile; then mv .Rprofile .Rprofile.DISABLED; fi
RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; fi
RUN make sass

#------------------------------------------------------------
# Expose port and set CMD
#------------------------------------------------------------
EXPOSE 3838
##CMD ["R", "-e", "shiny::runApp('components/app/R',port=3838,host='0.0.0.0')"]
CMD exec R -e "shiny::runApp('components/app/R', port=3838, host='0.0.0.0', launch.browser=0)" 2>&1 | tee -a run.log