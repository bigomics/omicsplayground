##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

multiwgcna_ai_report_controls_ui <- function(id) {
  multiomics_ai_report_controls_ui(id, module_label = "Module:")
}

multiwgcna_ai_report_controls_server <- function(id, module_choices = NULL) {
  multiomics_ai_report_controls_server(id, module_choices = module_choices)
}
