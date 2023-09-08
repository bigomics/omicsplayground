/*
##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##
*/


$(document).on('shiny:idle', function() {
    console.log("shiny:idle")
    Shiny.onInputChange("continue_test", true);
});

$(document).on('shiny:busy', function() {    
    Shiny.onInputChange("continue_test", false);
    console.log("shiny:busy")
});

