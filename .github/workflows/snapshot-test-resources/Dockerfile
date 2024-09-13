##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground:latest

#------------------------------------------------------------
# Set env variables 
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
ARG BRANCH_PLAYBASE=main
ARG GITHUB_SHA=master
ARG HEAD_REF=main
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen


# check if playbase@{$HEAD_REF} exists, if it does, install.. other wise install main
RUN echo "Checking playbase branch ${HEAD_REF}"
RUN BRANCH_EXISTS=$(git ls-remote --heads https://github.com/bigomics/playbase.git ${HEAD_REF} | wc -l) && \
    if [ $BRANCH_EXISTS -eq 1 ]; then \
    R -e "remotes::install_github('bigomics/playbase@${HEAD_REF}',dependencies=FALSE)"; \
    else \
    R -e "remotes::install_github('bigomics/playbase@${BRANCH_PLAYBASE}',dependencies=FALSE)"; \
    fi

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------

RUN rm -fr /omicsplayground  \
    && git clone -b ${BRANCH} --depth 1 \
    https://github.com/bigomics/omicsplayground.git /omicsplayground \
    && cd /omicsplayground \
    && git fetch origin $GITHUB_SHA \
    && git checkout $GITHUB_SHA \
    && chmod -R ugo+rwX /omicsplayground

WORKDIR /omicsplayground

RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; fi
RUN make sass

#------------------------------------------------------------
# Add testing framework dependencies
#------------------------------------------------------------
RUN R -e "install.packages('shinytest2')"

# add entrypoint to execute when the docker container starts up

# give permission to execute the entrypoint
RUN chmod +x .github/workflows/snapshot-test-resources/entrypoint.sh

ENTRYPOINT [".github/workflows/snapshot-test-resources/entrypoint.sh"]
