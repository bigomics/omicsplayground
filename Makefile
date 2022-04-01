
default:
	R -e "shiny::runApp(appDir='shiny',launch.browser=TRUE,port=3838)"

build.testing:
	docker build -f docker/Dockerfile.testing --no-cache -t bigomics/omicsplayground:testing . 

clean:
	rm `find -name '*~'`

run:
	docker run --rm -p 4000:3838 bigomics/omicsplayground:latest 

run.testing:
	docker run --rm -p 4000:3838 bigomics/omicsplayground:testing

run.develop:
	docker run --rm -p 4000:3838 bigomics/omicsplayground:develop

