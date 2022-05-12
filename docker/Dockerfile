##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

FROM bigomics/omicsplayground:base-ub2204
##FROM bigomics/omicsplayground:v2.3.0

#------------------------------------------------------------
# Install any extra (forgotten...) Ubuntu libs
#------------------------------------------------------------
ENV DEBIAN_FRONTEND noninteractive

RUN apt update && apt install -y locales

# Set the locale
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

#------------------------------------------------------------
# Install any extra (forgotten...) R packages since base
#------------------------------------------------------------

RUN R -e "remotes::install_github('bigomics/bigdash')"
# RUN R -e "install.packages('https://www.bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz')"
# RUN R -e "BiocManager::install(c('edgeR','DESeq2'))"
# RUN pip3 install umap-learn

#------------------------------------------------------------
# Download fresh code from GitHub
#------------------------------------------------------------
WORKDIR /
ARG TAG=master
RUN echo Downloading $TAG branch
RUN wget -nv https://github.com/bigomics/omicsplayground/archive/$TAG.zip \
    && rm -fr /omicsplayground*  \
    && unzip $TAG.zip \
    && mv omicsplayground-$TAG omicsplayground \
    && chmod -R ugo+rwX /omicsplayground \
    && rm $TAG.zip   

#------------------------------------------------------------
# Update configuration files into the Docker image
#------------------------------------------------------------
WORKDIR /omicsplayground
RUN mv .Rprofile .Rprofile.disable
RUN make sass

#------------------------------------------------------------
# Expose port and set CMD
#------------------------------------------------------------
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('components/app/R',port=3838,host='0.0.0.0')"]
