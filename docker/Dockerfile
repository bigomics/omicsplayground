##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground-base:ub2204_v3

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

RUN echo "Step 1: Install supporing libraries"

# This can be removed once new base is built
RUN R -e "remotes::install_github('cran/Matrix', ref = '2b05c7bd06c47eeb5939180938513da91d94b01b')"
RUN R -e "install.packages('plsRcox', dependencies = FALSE)"
# NEW INSTALLS THAT ARE BOTH IN OPG-BASE (can be removed when base v4 is built)

RUN echo "Step 2: Install supporing libraries"

# sf also added in requirements (check before removing!)
RUN R -e "install.packages('shinybrowser')"
RUN R -e "install.packages('sf')"

RUN echo "Step 3: Install supporting bigomics  (as R packages)"

# This can be removed once new base is built
RUN R -e "remotes::install_version('matrixStats',version='1.1.0',force=TRUE)"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz', type = 'source', repos = NULL)"
RUN R -e "BiocManager::install('orthogene', dependencies = NA, force = FALSE, ask = FALSE, update = FALSE)"

RUN echo "Step 4: Install supporing libraries"

RUN R -e "remotes::install_version('rgeos', version = '0.6-4')"
RUN R -e "remotes::install_version('SeuratObject', '4.1.4', repos = c('https://satijalab.r-universe.dev', getOption('repos')))"
RUN R -e "remotes::install_version('Seurat', '4.4.0', repos = c('https://satijalab.r-universe.dev', getOption('repos')))"
RUN R -e "BiocManager::install('org.Rn.eg.db', dependencies = NA, force = FALSE, ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('AnnotationHub',dependencies = NA, force = FALSE, ask = FALSE, update = FALSE)"
RUN R -e "devtools::install_version('dbplyr', version = '2.3.4')"

RUN echo "Step 5: Install supporing libraries"


WORKDIR /
##RUN R -e "remotes::install_github('bigomics/shinyChatR',dependencies=FALSE)"  ## rename??
RUN R -e "remotes::install_github('bigomics/playbase',dependencies=FALSE)"
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
