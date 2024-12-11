BRANCH:=`git rev-parse --abbrev-ref HEAD`  ## get active GIT branch
BRANCH:=$(strip $(BRANCH))
TAG=$(BRANCH)
ARG=

run: sass version rm.locks
	Rscript dev/run_app.R

run.headless:
	Rscript dev/run_app_headless.R

run.tee: 
	Rscript dev/run_app.R 2>&1 | tee -a run.log

run.profvis:
	Rscript dev/run_app_profvis.R

run.browser: rm.locks
	R --interactive < dev/run_app.R

sass: FORCE
	Rscript dev/sass.R
	Rscript dev/create_source_all.R

clean: rm.locks
	rm `find . -name '.#*' -o -name '*~' -o -name '#*#' -o -name '*.save' -o -name '*.bak' -o -name '*.SAVE' -o -name '*.BAK'`

rm.locks:
	find . -name 'LOCK*' -exec rm {} \;

show.branch:
	@echo $(BRANCH)

docker.run:
	@echo running docker *$(TAG)* at port 4000
	docker run --rm -it -p 4000:3838 bigomics/omicsplayground:$(TAG)

docker.run2:
	@echo running docker $(TAG) at port 4000
	docker run --rm -it -p 4000:3838 \
		-v ~/Playground/pgx:/omicsplayground/data \
		-v ~/Playground/libx-v4:/omicsplayground/libx \
		-v ~/Playground/pgx-visreport:/pgx-visreport \
		-v /aws/pgx-share:/omicsplayground/data_shared \
		-v /aws/pgx-public:/omicsplayground/data_public \
		-v ~/Playground/omicsplayground/etc:/omicsplayground/etc \
		bigomics/omicsplayground:$(TAG)

docker: FORCE version
	@echo building docker $(BRANCH)
	docker build $(ARG) --no-cache --build-arg BRANCH=$(BRANCH) \
		--progress plain \
		-f docker/Dockerfile \
	  	-t bigomics/omicsplayground:$(BRANCH) . \
		2>&1 | tee docker.log

update_playdata=false
update_bigdash=false
update_playbase=true
docker.update: FORCE
	@echo building update docker 
	docker build --no-cache \
		--progress plain \
		--build-arg update_playdata=$(update_playdata) \
		--build-arg update_bigdash=$(update_bigdash) \
		--build-arg update_playbase=$(update_playbase) \
		-f docker/Dockerfile.update \
	  	-t bigomics/omicsplayground:update .

docker.test: FORCE
	@echo building test docker 
	docker build --no-cache \
		--progress plain \
		-f docker/Dockerfile.test \
	  	-t bigomics/omicsplayground:test .

docker.bash:
	@echo bash into docker $(TAG)
	docker run -it bigomics/omicsplayground:$(TAG) /bin/bash

doc: FORCE
	Rscript dev/02_doc.R

install: FORCE
	Rscript dev/requirements.R

renv: FORCE
	R -e "renv::activate();renv::restore()"

FORCE: ;

DATE = `date +%y%m%d|sed 's/\ //g'`
VERSION = "v3.99.0"
BUILD := $(VERSION)"+"$(BRANCH)""$(DATE)

version: FORCE
	@echo "new version ->" $(BUILD)
	echo $(BUILD) > VERSION

changelog:
	sh ./dev/create-changelog.sh '.*' 3 >  CHANGELOG.md
	sh ./dev/create-changelog.sh '.*' 999 >  CHANGELOG-full.md
	sh ./dev/create-changelog.sh 'feat' 3 > FEATURES.md
	sh ./dev/create-changelog-pr.sh 1 4 > CHANGELOG-pr.md 

tags: version changelog
	git tag -f -a $(VERSION) -m 'version $(VERSION)'
	git push && git push --tags

push.latest: 
	docker tag bigomics/omicsplayground:$(BRANCH) bigomics/omicsplayground:latest
	docker push bigomics/omicsplayground:latest

push.version: 
	docker tag bigomics/omicsplayground:$(BRANCH) bigomics/omicsplayground:$(VERSION)
	docker push bigomics/omicsplayground:$(VERSION)

auth=none

board.launch:
	R -e "options(board = '$(board)', authentication = '$(auth)'); shiny::runApp('dev/board.launch')"

board.example:
	R -e "options(board = '$(board)', use_example_data = TRUE, authentication = '$(auth)'); shiny::runApp('dev/board.launch')"

pgx.check.error: sass
	Rscript dev/board_check_across_pgx.R $(if $(d),-d $(d),)

app.test: sass
	R -e  "options(authentication='$(auth)'); shiny::runTests()"

app.test.review:
	R -e "testthat::snapshot_review('snapshot/')"

update:
	Rscript dev/update_packages.R
