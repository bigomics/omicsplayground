##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

output <- sass::sass(
  sass::sass_file(
    'scss/main.scss'
  ),
  cache = NULL,
  options = sass::sass_options(
    output_style = "compressed"
  ),
  output = 'components/app/R/www/styles.min.css'
)
