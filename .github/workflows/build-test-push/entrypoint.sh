#!/bin/bash


echo "::set-output name=working_directory::$(pwd)"

# Run tests
R -e "x <- shiny::runTests(assert = FALSE); writeLines(as.character(all(x[[2]])), 'test_result.txt')"


# Read test results from file
test_result=$(cat test_result.txt)
echo "Test results: $test_result"

# return test result as an output
echo "test_result=$test_result" >> $GITHUB_OUTPUT