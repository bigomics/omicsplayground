
default:
	R -e "shiny::runApp(appDir='shiny',launch.browser=TRUE,port=3838)"

build.testing:
	docker build -f docker/Dockerfile.testing --no-cache -t bigomics/omicsplayground:testing . 

run.testing:
	docker run --rm -p 3838:3838 --name=playground bigomics/omicsplayground:testing

run:
	docker pull bigomics/omicsplayground:latest
	docker run --rm -p 3838:3838 --name=playground bigomics/omicsplayground:latest


