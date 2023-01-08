BRANCH=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch
TAG=latest
ifeq ($(BRANCH),'develop')
  TAG=develop
endif

VERSION=`head -n1 VERSION`

run:
	R -e "shiny::runApp('shiny',launch.browser=TRUE,port=3838)"

run.headless:
	R -e "shiny::runApp('shiny',launch.browser=FALSE,port=3838,host='0.0.0.0')"

clean:
	rm `find -name '*~'`

show.branch:
	@echo $(BRANCH)

run.docker:
	@echo running docker $(TAG) at port 4000
	docker run --rm -p 4000:3838 bigomics/omicsplayground:$(TAG)

run.docker2:
	@echo running docker $(TAG) at port 4000
	docker run --rm -p 4000:3838 -v /home/kwee/Playground/pgx:/omicsplayground/data -v /home/kwee/Playground/omicsplayground-dev/libx:/omicsplayground/libx bigomics/omicsplayground:$(TAG)

build.docker:
	@echo building docker $(TAG) from branch $(BRANCH) 
	docker build --no-cache --build-arg BRANCH=$(BRANCH) \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(TAG) .
build.base:
	@echo building ubuntu BASE docker 
	docker build --no-cache \
		-f docker/Dockerfile.base \
	  	-t bigomics/omicsplayground:base .
build.base2:
	@echo building ubuntu BASE docker 
	docker build \
		-f docker/Dockerfile.base \
	  	-t bigomics/omicsplayground:base .

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
	docker tag bigomics/omicsplayground:$(TAG) bigomics/omicsplayground:$(VERSION)

push.docker: 
	docker push bigomics/omicsplayground:$(VERSION)
	docker push bigomics/omicsplayground:latest
