#!/bin/bash

# Run tests
# run target sass in makefile

R -e "options(authentication='none'); x <- shiny::runTests(assert = FALSE); writeLines(as.character(all(x[[2]])), 'test_result.txt')"

# Read test results from file
res=$(cat test_result.txt)

# # return test result as an output (will be deprecated)
echo ::set-output name=test_result::$res

# line above should be switched to this, but for some reason it does not work
#echo "test_result=$res" >> $GITHUB_OUTPUT