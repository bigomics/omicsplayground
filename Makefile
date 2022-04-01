
default:
	R -e "shiny::runApp('shiny',launch.browser=TRUE,port=3838)"

clean:
	rm `find -name '*~'`

run:
	docker run --rm -p 4000:3838 bigomics/omicsplayground:latest 

run.testing:
	docker run --rm -p 4000:3838 bigomics/omicsplayground:testing

run.develop:
	docker run --rm -p 4000:3838 bigomics/omicsplayground:develop


build.testing:
	docker build -f docker/Dockerfile --no-cache -t bigomics/omicsplayground:testing . 

build.develop:
	docker build -f docker/Dockerfile --no-cache --build-arg TAG=develop -t bigomics/omicsplayground:develop . 
