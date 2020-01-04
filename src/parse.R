library(rmarkdown)

sources <- c("readme.Rmd", "MonteCarlo.Rmd")
src_path <- "E:/R/Learning/MCMC/src"
output_path <- "E:/R/Learning/MCMC"

for (s in sources) {
  render(file.path(src_path, s), output_dir = output_path,
         output_format = github_document(toc = TRUE, toc_depth = 2,
                                         html_preview = FALSE,
                                         pandoc_args = "--webtex"))
}