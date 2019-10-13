# Omics Playground: Explore Omics Data Freely

The Omics Playground is a user-friendly and interactive web-based platform 
for the analysis and visualization of transcriptomics and proteomics data. 
Currently the platform handles gene expression microarray, RNA-seq and 
LC-MS/MS proteomics data, and supports two species, human and mouse. The Omics 
Playground has been in particular devised to also support single cell RNA-seq 
data, as well as traditional gene expression experiments. 

The overview of the platform is shown in the figure below. It consists of
two main components. The first component addresses the data
importing and preprocessing, which includes preparing the input data, filtering,
normalising and precomputing of statistics for some analyses. The second part is
composed of the online interface, which supports the real-time visualisation and
interaction with users. The interface is subdivided into Basic and Expert modes
in order to provide a customisable experience suited to each user???s background.

.. figure:: figures/overview.png
    :align: center
    :width: 100%



# Installation

You can either run the platform from the source code, or download the docker image.


## Run from source code

Download the latest release of the source code (includes more data files) by cloning
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


## Run using the Docker file

The docker file of the platform is available on `Docker hub 
<https://www.docker.com/bigomics>`__.
Follow the steps below to set up a running platform from the docker file:

1. Pull the docker image using the command::

    docker pull bigomics/playground
2. Then run the docker with::

    docker run --rm -p 80:3838 bigomics/playground. 
3. Open ``localhost`` in your browser to run the platform.

.. note::

    The docker image requires about 5GB hard disk space.
    
    
