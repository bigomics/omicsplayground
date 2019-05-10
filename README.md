# Omics Playground

This is the source code for the Omics Playground, a visualization and analytics platform for multi-omics data. You can either clone this repository, download the latest release (includes more data files), or download our Docker image.

## Run from source code

1. Download or clone this repository. 
2. Be sure you have installed all necessary R packges by running the files `requirements.R` and `requirements2.R`.
3. In the `scripts` folder, run `run-all.R` to build the datasets. This can take a couple of hours.
4. Change into the `shiny` folder and run `R -e "rmarkdown::run()"`.

## Run using the Docker file

Pull the docker image (about 5Gb!) using the command `docker pull bigomics/playground`. Then run the docker with 
`docker run --rm -p 80:3838 bigomics/playground`. Then open `localhost` in your browser.

