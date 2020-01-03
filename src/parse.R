library(rmarkdown)

sources <- c("readme.Rmd", "MonteCarlo.Rmd")
src_path <- file.path(getwd(), "src")
output_path <- getwd()

for (s in sources) {
  render(file.path(src_path, s), output_dir = output_path,
         output_format = github_document(toc = TRUE, toc_depth = 2, html_preview = FALSE))
}