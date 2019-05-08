# Omics Playground

This is the source code for the Omics Playground, a visualization and analytics platform for multi-omics data. You can either clone this repository, download the latest release (includes more data files), or download our Docker image.

## From source code

1. Download or clone this repository. 
2. Be sure you you have installed all necessary R packges by running the file `requirements.R`
3. Change into the `shiny` folder and run `R -e "rmarkdown::run()"`

## Using the Docker file

Pull the docker image (warning: about 4Gb!) using the command `docker pull bigomics/playground`. Then run the docker with 
`docker run --rm -p 80:3838 playground`. Then open `localhost` in your browser.

