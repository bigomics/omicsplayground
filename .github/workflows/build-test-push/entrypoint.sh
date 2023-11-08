#!/bin/bash

# Run tests
R -e "x <- shiny::runTests(assert = FALSE); writeLines(as.character(all(x[[2]])), 'test_result.txt')"

# Read test results from file
test_result=$(cat test_result.txt)

# # return test result as an output (will be deprecated)
# echo ::set-output name=test_result::$test_result

# line above should be switched to this, but for some reason it does not work
echo "test_result=$time" >> $GITHUB_OUTPUT