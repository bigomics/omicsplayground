
default:
	cd shiny && R -e "shiny::runApp(launch.browser=TRUE,port=3939)"

testing:
	docker build -f docker/Dockerfile.testing --no-cache -t bigomics/omicsplayground:testing . 

clean:
	rm `find -name '*~'`
