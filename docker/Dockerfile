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
RUN apt update && apt install -y locales

# Set the locale (uncomment UTF8)
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

#------------------------------------------------------------
# Download fresh code from GitHub
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

#------------------------------------------------------------
# Install any extra (forgotten...) R packages since base
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN R -e "remotes::install_github('bigomics/bigdash')"
#RUN R -e "install.packages(c('RSpectra'))"
# RUN R -e "BiocManager::install(c('org.Rn.eg.db'))"
# RUN pip3 install umap-learn

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