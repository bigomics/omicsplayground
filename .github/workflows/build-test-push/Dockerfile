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
# Set env variables 
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen

#------------------------------------------------------------
# Install de-novo bigomics stuff (as R packages)
#------------------------------------------------------------

RUN R -e "remotes::install_github('bigomics/PCSF',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/playbase',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigdash',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/bigLoaders',dependencies=FALSE)"

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------

RUN rm -fr /omicsplayground  \
    && git clone -b $BRANCH --single-branch --depth 1 \
       https://github.com/mauromiguelm/omicsplayground.git \   
    && chmod -R ugo+rwX /omicsplayground 

WORKDIR /omicsplayground

RUN if test -f renv.lock; then mv renv.lock renv.lock.DISABLED; fi
RUN make sass


# add entrypoint to execute when the docker container starts up
#RUN cp /.github/workflows/build-test-push/entrypoint.sh .
#RUN chmod +x /entrypoint.sh
#RUN cat /entrypoint.sh

# give permission to execute the entrypoint
RUN chmod +x .github/workflows/build-test-push/entrypoint.sh

ENTRYPOINT [".github/workflows/build-test-push/entrypoint.sh"]
