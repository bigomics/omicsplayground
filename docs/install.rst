.. _Installation:

Installation
================================================================================

The current version of the Omics Playground software is implemented in R 
using the `Shiny <https://shiny.rstudio.com/>`__ web application framework. 
You can either run the platform from the source code, or download the docker image.


Run from source code
--------------------------------------------------------------------------------
The source code of the platform is available on 
`GitHub <https://github.com/IRB-Bioinformatics/OmicsPlayground>`__. You can 
download the latest release of the software (includes more data files) by cloning
the repository. Below, we explain the steps required to set up the platform from
the source code:

1. Download or clone the GitHub repository to a location in your ``PATH``::

    git clone https://github.com/bigomics/playground.git
2. Be sure you have installed all necessary R packages by running the files in the ``/R`` folder::

    R requirements.R
    R requirements2.R
3. Similarly, run the following command in the ``/scripts`` folder to build the datasets::

    R run-all.R

.. note::

    Building the datasets can vary from minutes to a couple of hours depending on their sizes.
4. Change the current directory into the ``/shiny`` folder and execute the following command to run the platform::

    R -e "rmarkdown::run()"


Run using the Docker file
--------------------------------------------------------------------------------
The docker file of the platform is available on `Docker hub 
<https://cloud.docker.com/u/bigomics/repository/docker/bigomics/playground>`__.
Follow the steps below to set up a running platform from the docker file:

1. Pull the docker image using the command::

    docker pull bigomics/playground
2. Then run the docker with::

    docker run --rm -p 80:3838 bigomics/playground. 
3. Open ``localhost`` in your browser to run the platform.

.. note::

    The docker image requires about 5GB hard disk space.
    
    
