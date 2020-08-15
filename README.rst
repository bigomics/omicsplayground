Omics Playground: Explore Omics Data Freely
================================================================================

Omics Playground is a comprehensive *self-service analytics platform* for
the visualization, analytics and exploration of Big Omics Data. It allows
biologists to apply a multitude of state-of-the-art analysis tools to their
own data to explore and discover underlying biology without coding.

Installation
================================================================================

You can either run the platform from the source code, or download the
docker image. Running Omics Playground from the docker file is the
easiest way.
    
Run using the Docker file
--------------------------------------------------------------------------------
The docker file of the platform is available on `Docker Hub 
<https://hub.docker.com/r/bigomics/omicsplayground>`__.
Follow the steps below to set up a running platform from the docker file:

1. Pull the docker image using the command::

    docker pull bigomics/omicsplayground
    
   Warning. The docker image requires about 5GB-8GB hard disk space. Note: download
   version v1.0 if you want the exact version of the NAR/GAB publication, otherwise
   docker will download the latest version by default.
    
2. Now run the docker with::

    docker run --rm -p 4000:3838 bigomics/omicsplayground
    
3. Then open ``http://localhost:4000`` in your browser to access the platform.

   
   
Run from source code
--------------------------------------------------------------------------------

Omics Playground is implemented in R using the `Shiny Server
<https://shiny.rstudio.com/>`__ web application framework. You will
need R and Shiny Server installed to run Omics Playground. The source code of the platform is available on 
`GitHub <https://github.com/bigomics/omicsplayground>`__. You can 
download the latest release of the software by cloning the repository. 

Below, we explain the steps required to set up the platform from
the source code:

1. Clone the GitHub repository using::

    git clone https://github.com/bigomics/omicsplayground.git
   
   Note: download version v1.0 if you want the exact version of the NAR/GAB publication, 
   otherwise GitHub will download the latest version by default.
    
2. Install all necessary R packages by running the script in the ``R/`` folder::

    Rscript requirements.R
    
3. Run the following command in the ``build/`` folder to build the datasets::

    Rscript build-datasets.R

   Note: Building the datasets can vary from minutes to a couple of hours depending on their sizes.

4. Change the current directory into the ``shiny/`` folder and execute the following command
   to run the platform::

    R -e "shiny::runApp(launch.browser=TRUE)"

   If you have Shiny Server installed you can create a link to the
   shiny folder in the system-wide shiny-server apps folder or in your
   ShinyApps user folder.



Documentation
=======================================================================================

The platform consists of two main components. The first component is off-line and addresses the data
importing and preprocessing, which includes preparing the input data, filtering, 
normalising and precomputing of statistics for some analyses. The second part is
composed of the online interface, which supports the real-time visualisation and
interaction with users. The interface is subdivided into Basic and Expert modes
to provide a customisable experience suited to each user's background.

The docker image and the installation script will contain some example data sets. To analyze your
own data you can use the upload function, or create/modify the scripts in the ``scripts/``folder.
Creating a custom script is much more flexible and allows, if necessary, batch correction, 
quality filtering and/or translation of probe names.

More detailed information and feature explanation of Omics Playground is 
available in the `online documentation <https://omicsplayground.readthedocs.io>`__.

.. figure:: docs/figures/overview.png
    :align: center
    :width: 90%

