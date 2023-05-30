BRANCH=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch
BRANCH:=$(strip $(BRANCH))

run: sass tag.version
	R -e "shiny::runApp('components/app/R',launch=TRUE,port=3838)"

run.headless:
	R -e "shiny::runApp('components/app/R',launch=FALSE,port=3838,host='0.0.0.0')"

run.tee: 
	R -e "shiny::runApp('components/app/R',launch=TRUE,port=3838)" 2>&1 | tee -a run.log

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


TAG=$(BRANCH)
docker.run:
	@echo running docker *$(TAG)* at port 4000
	docker run --rm -it -p 4000:3838 bigomics/omicsplayground:$(TAG)

docker.run2:
	@echo running docker $(TAG) at port 4000
	docker run --rm -it -p 4000:3838 \
		-v ~/Playground/pgx:/omicsplayground/data \
		-v ~/Playground/libx:/omicsplayground/libx \
		-v ~/Playground/config/firebase:/omicsplayground/components/app/R/firebase \
		bigomics/omicsplayground:$(TAG)

docker: FORCE tag.version
	@echo building docker $(BRANCH)
	docker build --no-cache --build-arg BRANCH=$(BRANCH) \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(BRANCH) .

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
	@echo bash into docker $(BRANCH)
	docker run -it -p 3838:3838 bigomics/omicsplayground:$(BRANCH) /bin/bash


doc: FORCE
	Rscript dev/02_doc.R

install: FORCE
	Rscript dev/03_install.R

renv: FORCE
	R -e "renv::activate();renv::restore()"

FORCE: ;

##VERSION=`head -n1 VERSION`
DATE = `date +%y%m%d|sed 's/ //g'`
VERSION := "v3.0-RC11."$(BRANCH)""$(DATE)

tag.version:
	@echo "new version ->" $(VERSION)
	sed -i "1s/.*/$(VERSION)/" VERSION

tags:
	git tag -f -a $(VERSION) -m 'version $(VERSION)'
	git push && git push --tags

push.latest: 
	docker tag bigomics/omicsplayground:$(BRANCH) bigomics/omicsplayground:latest
	docker push bigomics/omicsplayground:latest

push.version: 
	docker tag bigomics/omicsplayground:$(BRANCH) bigomics/omicsplayground:$(VERSION)
	docker push bigomics/omicsplayground:$(VERSION)

BOARD = dataview
board:
	R -e "source('components/golem_utils/run_dev.R');launch_board('board.$(BOARD)', options=list(launch.browser=TRUE))"

