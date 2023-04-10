TAG=latest-v2
BRANCH=v2
VERSION=`head -n1 VERSION`

run:
	R -e "shiny::runApp('shiny',launch.browser=TRUE,port=3939)"

run.headless:
	R -e "shiny::runApp('shiny',launch.browser=FALSE,port=3939,host='0.0.0.0')"

clean:
	rm `find -name '*~'`

run.docker:
	@echo running docker $(TAG) at port 4000
	docker run --rm -it -p 4000:3838 bigomics/omicsplayground:$(TAG)

run.docker2:
	@echo running docker $(TAG) at port 4000
	docker run --rm -it -p 4000:3838 -v /home/kwee/Playground/pgx:/omicsplayground/data -v /home/kwee/Playground/omicsplayground-dev/libx:/omicsplayground/libx bigomics/omicsplayground:$(TAG)

run.docker3:
	@echo running docker $(TAG) at port 4000
	docker run --rm -it -p 4000:3838 -v /home/kwee/Playground/config/firebase.rds:/omicsplayground/shiny/firebase.rds -v /home/kwee/Playground/config/Renviron.site:/etc/R/Renviron.site -v /home/kwee/Playground/config/OPTIONS.fb:/omicsplayground/shiny/OPTIONS bigomics/omicsplayground:$(TAG)

build.docker:
	@echo building docker from branch $(BRANCH)
	docker build --no-cache \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(TAG) .

build.base:
	@echo building ubuntu BASE docker 
	docker build --no-cache \
		-f docker/Dockerfile.base \
	  	-t bigomics/omicsplayground-base:21.04 .

build.test:
	@echo building test docker 
	docker build --no-cache \
		-f docker/Dockerfile.test \
	  	-t bigomics/omicsplayground:test .

bash.docker:
	@echo bash into docker $(TAG)
	docker run -it bigomics/omicsplayground:$(TAG) /bin/bash

tags:
	git tag -f -a $(VERSION) -m 'version $(VERSION)'
	git push && git push --tags
	docker tag bigomics/omicsplayground:latest bigomics/omicsplayground:$(VERSION)

push.docker: 
	docker push bigomics/omicsplayground:$(VERSION)
	docker push bigomics/omicsplayground:latest
