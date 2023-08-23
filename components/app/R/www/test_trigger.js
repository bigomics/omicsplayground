/*
##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##
*/

// code included inside $(document).ready() will only run once the page is ready for JavaScript code to execute
$(document).ready(function() {
  
  // initialize a counter
  var continue_test = false;
  
  $("#logo-bigomics").on("click", function(){
  
    continue_test = true;
    
    // send message to Shiny
    Shiny.onInputChange("continue_test", continue_test);
  });

});