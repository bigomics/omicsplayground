#!/bin/bash

# Run Docker container 
docker run --rm -it -p 4000:3838 --name test_container bigomics/omicsplayground:latest

# capture test results
docker exec -it test_container R -e "shiny::runTests()" > test_result.txt

# Read test results from file
test_result=$(cat test_result.txt)
echo "Test results: $test_result"

# Check test results and run command
if [[ $test_result == FALSE ]]; then
  echo "All tests passed. Running command..."
  exec "$@"
else
  echo "Some tests failed. Not running command."
  exit 1
fi