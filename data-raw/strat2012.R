## code to prepare `strat2012` dataset

devtools::install_github("fmichonneau/phyloch")
utils::data("strat2012", package = "phyloch")
usethis::use_data(strat2012, overwrite = TRUE)
