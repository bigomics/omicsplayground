##
## compute-comparer: prove the Shiny upload flow and the edgy script path
## produce the same pgx. See README.md.
##
SHELL := /bin/bash
HERE  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

OPG_MAIN ?= $(abspath $(HERE)/../..)
OPG_EDGY ?= /home/massagno/bigomics/GitHub/omicsplayground-edgy
PB_MAIN  := $(HERE)/libs/pb-main
PB_EDGY  := $(HERE)/libs/pb-edgy

SC   ?= baseline
PORT ?= 3838
SEED ?= 42
RUN  := $(HERE)/runs/$(SC)

# R_LIBS_USER REPLACES the personal library rather than adding to it, so setting it
# bare hides every package installed only in ~/R/<arch>-library (imputeLCMD, needed
# by QRILC, is one -- it fails with "requireNamespace(\"imputeLCMD\") is not TRUE").
# Append the personal library so the harness playbase still wins on precedence but
# nothing else disappears.
R_USER_LIB := $(shell Rscript -e 'cat(Sys.getenv("R_LIBS_USER"))' 2>/dev/null)
R_PERSONAL := $(HOME)/R/x86_64-pc-linux-gnu-library/4.3

A_ENV := R_LIBS_USER=$(PB_MAIN):$(R_PERSONAL)
B_ENV := R_LIBS_USER=$(PB_EDGY):$(R_PERSONAL)

.PHONY: setup check app.start app.stop a b create.a create.b full.a full.b \
        g1 g2 g3 gate1 gate2 clean.runs clean.userinput

## ---------------------------------------------------------------- setup
setup:
	@mkdir -p $(HERE)/libs/pb-main $(HERE)/libs/pb-edgy $(HERE)/runs
	R CMD INSTALL -l $(HERE)/libs/pb-main /home/massagno/bigomics/GitHub/playbase-main
	R CMD INSTALL -l $(HERE)/libs/pb-edgy /home/massagno/bigomics/GitHub/playbase
	# BOTH playbase versions call playdata::PHOSPHOSITE, but the system playdata
	# (0.1.0, branch feat/reactome-diagrams) does not export it -- createPGX dies
	# with "'PHOSPHOSITE' is not an exported object". playdata@main has it. Install
	# into the harness libs only, so the user's system playdata stays as-is.
	R CMD INSTALL -l $(HERE)/libs/pb-main /home/massagno/bigomics/GitHub/playdata-main
	R CMD INSTALL -l $(HERE)/libs/pb-edgy /home/massagno/bigomics/GitHub/playdata-main
	cd $(HERE) && npm i playwright

check:
	@echo -n "pb-main : "; $(A_ENV) Rscript -e 'suppressMessages(library(playbase)); cat(as.character(packageVersion("playbase")), "preprocess=", "pgx.preprocess" %in% getNamespaceExports("playbase"), "\n")' 2>/dev/null | tail -1
	@echo -n "pb-edgy : "; $(B_ENV) Rscript -e 'suppressMessages(library(playbase)); cat(as.character(packageVersion("playbase")), "preprocess=", "pgx.preprocess" %in% getNamespaceExports("playbase"), "\n")' 2>/dev/null | tail -1
	@echo "opg main: $(OPG_MAIN)"
	@echo "opg edgy: $(OPG_EDGY)"

## ---------------------------------------------------------------- the app
## Path A's app MUST run against playbase main out of the master worktree.
app.start:
	@cd $(OPG_MAIN) && $(A_ENV) nohup Rscript dev/run_app_headless.R \
		> $(HERE)/runs/app.log 2>&1 & echo "app starting, log: runs/app.log"
	@until curl -sf -o /dev/null http://localhost:$(PORT); do sleep 2; done; echo "app up on :$(PORT)"

app.stop:
	@pkill -f "run_app_headless.R" 2>/dev/null; echo "app stopped"

## ------------------------------------------------------------- capture A/B
## A: drive the browser, capture params.RData before the app deletes raw_dir.
a:
	cd $(HERE) && node drive_app.mjs $(SC) 2>&1 | tee runs/driver-$(SC).log

## B: build the script-side params from the same fixture + scenario.
b:
	@mkdir -p $(RUN)/B
	cd $(HERE) && $(B_ENV) Rscript run_script.R $(SC) 2>&1 | tee runs/script-$(SC).log

## ---------------------------------------------------------------- execute
create.a:
	$(A_ENV) Rscript $(HERE)/run_pgx.R $(RUN)/A $(OPG_MAIN) create $(SEED)
create.b:
	$(B_ENV) Rscript $(HERE)/run_pgx.R $(RUN)/B $(OPG_EDGY) create $(SEED)

full.a:
	$(A_ENV) Rscript $(HERE)/run_pgx.R $(RUN)/A $(OPG_MAIN) full $(SEED)
full.b:
	$(B_ENV) Rscript $(HERE)/run_pgx.R $(RUN)/B $(OPG_EDGY) full $(SEED)

## ------------------------------------------------------------------ gates
## Gate 1 needs pgx.preprocess, so it runs under pb-edgy.
g1:
	$(B_ENV) Rscript $(HERE)/compare.R $(SC) 1
g2:
	$(B_ENV) Rscript $(HERE)/compare.R $(SC) 2
g3:
	$(B_ENV) Rscript $(HERE)/compare.R $(SC) 3

## Convenience: the two cheap gates end to end for one scenario.
gate1: a b g1
gate2: create.a create.b g2

## ----------------------------------------------------------------- hygiene
## The poller keys on raw_ dirs that appear after compute is triggered, but a
## stale one from a crashed run is still confusing. Clear between batches.
clean.userinput:
	rm -rf $(OPG_MAIN)/data/USER_INPUT/raw_*

clean.runs:
	rm -rf $(HERE)/runs/*/ $(HERE)/runs/report.csv
