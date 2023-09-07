#!/bin/bash

# Run tests
R --slave -e "cat('Current working directory:', getwd(), '\n')"
R -e "shiny::runTests();all(shiny::runTests()[[2]])"
R --slave -e "writeLines(all(shiny::runTests()[[2]]), 'test_result.txt')"
#R --slave -e "all(shiny::runTests()[[2]])" > test_result.txt

# Read test results from file
test_result=$(cat test_result.txt)
echo "Test results: $test_result"

# return test result as an output
echo ::set-output name=test_result::$test_result

# Check test results and run command
if [[ $test_result == FALSE ]]; then
  echo "All tests passed. Running command..."
  exec "$@"
else
  echo "Some tests failed. Not running command."
  exit 0
fi