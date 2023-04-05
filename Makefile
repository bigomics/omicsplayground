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

run.tee: 
	R -e "shiny::runApp('components/app/R',launch=TRUE,port=3838)" 2&>1 | tee -a run.log

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

docker.run:
	@echo running docker $(TAG) at port 4000
	docker run --rm -it -p 4000:3838 bigomics/omicsplayground:$(TAG)

docker.run2:
	@echo running docker $(TAG) at port 4001
	docker run --rm -it -p 4001:3838 bigomics/omicsplayground:$(TAG)

docker: FORCE
	@echo building docker $(TAG)
	docker build --no-cache --build-arg TAG=$(TAG) \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(TAG) .

docker.base: FORCE
	@echo building docker BASE
	docker build  \
		-f docker/Dockerfile.base \
	  	-t bigomics/omicsplayground-base:ub2204 .

docker.test: FORCE
	@echo building test docker 
	docker build --no-cache \
		-f docker/Dockerfile.test \
	  	-t bigomics/omicsplayground:test .

docker.bash:
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

date = `date +%y%m%d|sed 's/ //g'`
version = "v3.0-rc3".$(date).$(BRANCH)
version:
	@echo "new version ->" $(version)
#	sed -i "1s/.*/$(version)/" VERSION
