
Omics Playground: Explore Omics Data Freely
================================================================================

Omics Playground is a comprehensive *self-service analytics platform* for
the visualization, analytics and exploration of Big Omics Data. It allows
biologists to apply a multitude of state-of-the-art analysis tools to their
own data to explore and discover underlying biology without coding.

How to use Omics Playground?
=======================================================================================
The detailed documentation on how to load data and how to use the functionalities of omics playground can be found `here <https://omicsplayground.readthedocs.io>`__ for text-based tutorials, or `here <https://bigomics.ch/tutorials/>`__ for video tutorials.


Platform components
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

Run from source code / Start Developing
--------------------------------------------------------------------------------

Omics Playground relies on 3 basic components that do the work on the background: playbase, bigdash, bigloaders. Thus, it is necessary install them manually within the ``R`` environmnet::

    remotes::install_github('bigomics/playbase')
    remotes::install_github('bigomics/bigdash')
    remotes::install_github('bigomics/bigLoaders')

On top of these, a python interpreter is also necessary for the interactive plots. This can be aslo easily installed all within R via::

    install.packages("reticulate")
    reticulate::install_miniconda()

Then, everything is ready for installing omicsplayground::

    git clone https://github.com/bigomics/omicsplayground.git
   
Note: download version v1.0 if you want the exact version of the NAR/GAB publication, otherwise GitHub will download the latest version by default.
    
Next, install all necessary R packages and dependencies by running from the omicsplayground folder::

    cd omicsplayground
    Rscript dev/requirements.R
    
Finally, you can run the omicsplayground platform. You can do this with the Makefile located in the root omicsplayground folder::

    make run
    
Or you can launch the platform from within an R session::

   shiny::runApp('components/app/R', launch.browser=TRUE)

   If you have Shiny Server installed you can create a link to the
   shiny folder in the system-wide shiny-server apps folder or in your
   ShinyApps user folder.
