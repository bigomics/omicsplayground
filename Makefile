BRANCH=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch
TAG=latest
ifeq ($(BRANCH),'develop')
  TAG=develop
endif
VERSION=`head -n1 VERSION`

run: sass
	R -e "shiny::runApp('components/app/R',launch=TRUE,port=3838)"

run.headless:
	R -e "shiny::runApp('components/app/R',launch=FALSE,port=3838,host='0.0.0.0')"

sass: FORCE
	Rscript dev/sass.R
	Rscript dev/create_source_all.R

clean:
	find . -name '.#*' -o -name '#*' -o -name '*~'
	@echo -n "Are you sure? [y/N] " && read ans && [ $${ans:-N} = y ]
	rm `find . -name '.#*' -o -name '#*' -o -name '*~'`

clean.force:
	rm `find . -name '.#*' -o -name '#*' -o -name '*~'`

show.branch:
	@echo $(BRANCH)

run.docker:
	@echo running docker $(TAG) at port 3838
	docker run --rm -p 3838:3838 bigomics/omicsplayground:$(TAG)

build.base:
	@echo building docker BASE
	docker build --no-cache \
		-f docker/Dockerfile.base \
	  	-t bigomics/omicsplayground:base .

build.docker:
	@echo building docker $(TAG)
	docker build --no-cache --build-arg TAG=$(TAG) \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(TAG) .

bash.docker:
	@echo bash into docker $(TAG)
	docker run -it -p 3838:3838 bigomics/omicsplayground:$(TAG) /bin/bash


doc: FORCE
	Rscript dev/02_doc.R

install: FORCE
	Rscript dev/03_install.R

renv: FORCE
	R -e "renv::activate();renv::restore()"

FORCE: ;

tags:
	git tag -f -a $(VERSION) -m 'version $(VERSION)'
	git push && git push --tags

push.latest: 
	docker tag bigomics/omicsplayground:testing bigomics/omicsplayground:latest
	docker tag bigomics/omicsplayground:testing bigomics/omicsplayground:$(VERSION)
	docker push bigomics/omicsplayground:latest
	docker push bigomics/omicsplayground:$(VERSION)
