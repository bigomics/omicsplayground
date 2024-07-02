##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground-base:ub2204_v4

ARG BRANCH=master
RUN echo "Building Docker image from branch $BRANCH"

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
COPY docker/Rprofile .Rprofile
#RUN pip3 install umap-learn

#------------------------------------------------------------
# Install supporing libraries
#------------------------------------------------------------
RUN R -e "remotes::install_version('Matrix', version = '1.6.1.1')"
RUN R -e "install.packages('ruv')"


WORKDIR /

# check if playbase@{$BRANCH} exists, if it does, install.. other wise install main
RUN echo "Checking playbase branch ${BRANCH}"
RUN BRANCH_EXISTS=$(git ls-remote --heads https://github.com/bigomics/playbase.git ${BRANCH} | wc -l) && \
    if [ $BRANCH_EXISTS -eq 1 ]; then \
    R -e "remotes::install_github('bigomics/playbase@${BRANCH}',dependencies=FALSE)"; \
    else \
    R -e "remotes::install_github('bigomics/playbase@main',dependencies=FALSE)"; \
    fi
RUN R -e "remotes::install_github('bigomics/bigdash',dependencies=FALSE)"
RUN R -e "remotes::install_github('the-y-company/bsutils',dependencies=FALSE, force=TRUE)"
RUN R -e "remotes::install_github('bigomics/wizardR',dependencies=FALSE, force=TRUE)"


#------------------------------------------------------------
# Download fresh code from GitHub (not an R package)
#------------------------------------------------------------
WORKDIR /
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
