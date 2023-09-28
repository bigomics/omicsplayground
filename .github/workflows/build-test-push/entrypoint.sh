#!/bin/bash

# Set output variable with current working directory
echo "::set-output name=working_directory::$(pwd)"

# Run tests
R -e "x <- shiny::runTests(assert = FALSE); writeLines(as.character(all(x[[2]])), 'test_result.txt')"

# Read test results from file
test_result=$(cat test_result.txt)

# Set output variable with test results
echo "::set-output name=test_result::$test_result"

# Write test result to file
echo "$test_result" > test_result_output.txt

# Return test result as an output
echo "test_result=$test_result" >> $GITHUB_OUTPUT