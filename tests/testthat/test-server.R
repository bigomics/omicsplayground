
# Configure this test to fit your need.
# testServer() function makes it possible to test code in server functions and modules, without needing to run the full Shiny application
#testServer(app_server, {

  # Set and test an input
#  session$setInputs(x = 2)
#  expect_equal(input$x, 2)

  # Example of tests you can do on the server:
  # - Checking reactiveValues
  # expect_equal(r$lg, 'EN')
  # - Checking output
  # expect_equal(output$txt, "Text")
#})

# Configure this test to fit your need
#test_that(
#  "app launches",
#  {
#    golem::expect_running(sleep = 5)
#  }
#)
