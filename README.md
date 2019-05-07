# Omics Playground

This is the source code for the Omics Playground, an R/Shiny analysis and visualization platform for multi-omics data. You can either run the app from the prepared Docker file, or by cloning this repository.


## Using the Docker file

Pull the docker image (warning: it's about 4Gb) using the command `docker pull bigomics/playground`. Then run the docker as 
`docker run --rm -p 80:3838 playground`. Then open `localhost` in your browser.


## From source code

1. Download or clone this repository. 
2. Be sure you you have installed all necessary R packges by running the file `require.R`
3. Change into the `shiny` folder and run `R -e "rmarkdown::run()"`
