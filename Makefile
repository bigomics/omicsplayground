BRANCH=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch
TAG=latest
ifeq ($(BRANCH),'develop')
  TAG=develop
endif

run:
	R -e "shiny::runApp('components/app/R',launch.browser=TRUE,port=3838)"

run.headless:
	R -e "shiny::runApp('shiny',launch.browser=FALSE,port=3838,host='0.0.0.0')"

clean:
	find . -name '.#*' -o -name '#*' -o -name '*~'
	@echo -n "Are you sure? [y/N] " && read ans && [ $${ans:-N} = y ]
	rm `find . -name '.#*' -o -name '#*' -o -name '*~'`

clean.force:
	rm `find . -name '.#*' -o -name '#*' -o -name '*~'`

show.branch:
	@echo $(BRANCH)

run.docker:
	@echo running docker $(TAG) at port 4000
	docker run --rm -p 4000:3838 bigomics/omicsplayground:$(TAG)

build.base:
	@echo building docker BASE
	docker build --no-cache \
		-f dev/docker/Dockerfile.base \
	  	-t bigomics/omicsplayground:base .

build.docker:
	@echo building docker $(TAG)
	docker build --no-cache --build-arg TAG=$(TAG) \
		-f dev/docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(TAG) .

bash.docker:
	@echo bash into docker $(TAG)
	docker run -it bigomics/omicsplayground:$(TAG) /bin/bash


doc: FORCE
	Rscript dev/02_doc.R

install: FORCE
	Rscript dev/03_install.R


FORCE: ;
