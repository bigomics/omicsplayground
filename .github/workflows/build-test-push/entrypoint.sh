#!/bin/bash

# Run tests
R -e "shiny::runTests(assert = FALSE);all(shiny::runTests(assert = FALSE)[[2]])"
R --slave -e "writeLines(as.character(all(shiny::runTests(assert = FALSE)[[2]])), 'test_result.txt')"

# Read test results from file
test_result=$(cat test_result.txt)
echo "Test results: $test_result"

# return test result as an output
echo ::set-output name=test_result::$test_result
echo "{test_result}={$test_result}" >> $test_result