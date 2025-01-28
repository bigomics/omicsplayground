##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from playbase image
#------------------------------------------------------------

FROM bigomics/playbase:latest

ARG BRANCH=master
RUN echo Building Docker image from branch '$BRANCH'

ENV DEBIAN_FRONTEND=noninteractive
# Set the locale (uncomment UTF8)
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

#------------------------------------------------------------
# Install any extra Ubuntu libs (should be moved to playbase)
#------------------------------------------------------------
RUN apt update && apt install -y fonts-lato

WORKDIR /tmp

RUN ARCH=$(dpkg --print-architecture) && \
    if [ "$ARCH" = "arm64" ]; then \
        wget -nv https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.56/quarto-1.5.56-linux-arm64.deb; \
    elif [ "$ARCH" = "amd64" ]; then \
        wget -nv https://github.com/quarto-dev/quarto-cli/releases/download/v1.5.56/quarto-1.5.56-linux-amd64.deb; \
    else \
        echo "Unsupported architecture: $ARCH"; exit 1; \
    fi && \
    dpkg -i quarto-1.5.56-linux-*.deb; \
    quarto install tinytex

#------------------------------------------------------------
# Add supporting R libraries
#------------------------------------------------------------
WORKDIR /
COPY dev/Rprofil[e] /root/.Rprofile
# RUN R -e "install.packages('plotly')"

#------------------------------------------------------------
# Download fresh omicsplayground code from GitHub
#------------------------------------------------------------

RUN echo Downloading $BRANCH branch
RUN rm -fr /omicsplayground  \
    && git clone -b $BRANCH --single-branch --depth 1 \
    https://github.com/bigomics/omicsplayground.git \   
    && chmod -R ugo+rwX /omicsplayground 

# Run auto-install script
WORKDIR /omicsplayground
COPY dev/*.R dev/description* dev/
RUN Rscript dev/requirements.R

#------------------------------------------------------------
# Install fresh playbase (from same BRANCH)
#------------------------------------------------------------
# check if playbase@{$BRANCH} exists, if it does, install.. otherwise install main
RUN echo "Checking playbase branch ${BRANCH}"
RUN BRANCH_EXISTS=$(git ls-remote --heads https://github.com/bigomics/playbase.git ${BRANCH} | wc -l) && \
    if [ $BRANCH_EXISTS -eq 1 ]; then \
        R -e "remotes::install_github('bigomics/playbase@${BRANCH}',dependencies=FALSE)"; \
    else \
        R -e "remotes::install_github('bigomics/playbase@main',dependencies=FALSE)"; \
    fi

# after install non-detected missing packages
#RUN R -e "remotes::install_version('Matrix', version = '1.6.1.1')"
RUN R -e "install.packages(c('shiny.i18n','waiter','sever'))"
#RUN R -e "remotes::install_github('bigomics/bigdash',dependencies=FALSE)"
#RUN R -e "remotes::install_github('the-y-company/bsutils',dependencies=FALSE)"

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------

## initialize
WORKDIR /omicsplayground
RUN make sass
#RUN Rscript dev/init_annothub.R

## disable Rprofile and Renv 
RUN if test -f .Rprofile; then mv .Rprofile .Rprofile.DISABLED; fi
RUN if test -f .Renviron; then mv .Renviron .Renviron.DISABLED; fi
RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; fi

#------------------------------------------------------------
# Expose port and set entry CMD
#------------------------------------------------------------
EXPOSE 3838
CMD exec Rscript dev/run_app_headless.R 2>&1 | tee -a run.log
