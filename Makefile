TAG=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch

run:
	R -e "shiny::runApp('shiny',launch.browser=TRUE,port=3838)"

run.headless:
	R -e "shiny::runApp('shiny',launch.browser=FALSE,port=3838,host='0.0.0.0')"

clean:
	rm `find -name '*~'`

show.tag:
	@echo $(TAG)

run.docker:
	@echo running docker $(TAG) at port 4000
	docker run --rm -p 4000:3838 bigomics/omicsplayground:$(TAG)

build.docker:
	@echo building docker $(TAG)
	docker build --no-cache --build-arg TAG=$(TAG) \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(TAG) .

bash.docker:
	@echo bash into docker $(TAG)
	docker run -it bigomics/omicsplayground:$(TAG) /bin/bash

