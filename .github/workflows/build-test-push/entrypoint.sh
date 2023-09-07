#!/bin/bash

# Run tests
R -e "shiny::runTests();all(shiny::runTests()[[2]])"
R --slave -e "writeLines(as.character(all(shiny::runTests()[[2]])), 'test_result.txt')"
#R --slave -e "all(shiny::runTests()[[2]])" > test_result.txt

# Read test results from file
test_result=$(cat test_result.txt)
echo "Test results: $test_result"

# return test result as an output
echo ::set-output name=test_result::$test_result