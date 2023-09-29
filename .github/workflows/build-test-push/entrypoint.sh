#!/bin/bash

echo "::set-output name=working_directory::$(sudo pwd)"
# Run tests
R -e "x <- shiny::runTests(assert = FALSE); writeLines(as.character(all(x[[2]])), 'test_result.txt')"

# Read test results from file
test_result=$(cat test_result.txt)
echo "Test results: $test_result"

# return git diff to output
git_diff=$(git diff)
echo "::set-output name=git_diff::$(git_diff)"

# return test result as an output
echo ::set-output name=test_result::$test_result
#echo "{test_result}={$test_result}" >> $test_result # return test result as an output

