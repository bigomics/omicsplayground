R Style Guide
================

R is a high-level programming language used primarily for statistical
computing and graphics. The goal of the R Programming Style Guide is
to make our R code easier to read, share, and verify. Our style guide is
largely based on the following R Style Guides:

1. [Google's R Style Guide](http://web.stanford.edu/class/cs109l/unrestricted/resources/google-style.html)
2. [The tidyverse style guide](https://style.tidyverse.org/)
3. [Hadley Wickahm's style guid](http://adv-r.had.co.nz/Style.html)


File Names
---------------------

File names should end in .R and be meaningful. Some older files still
end with *.r. Files in the folder R are prefixed by category (gx-,
ngs-, pgx-) that mostly correspond to the analysis type and base input
object used.

- GOOD: predict_ad_revenue.R
- BAD: foo.R

Identifiers
-----------------------
