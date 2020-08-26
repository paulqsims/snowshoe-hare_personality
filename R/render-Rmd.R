################################################################################
# Purpose: Markdown and HTML generation for Rmd files for Lafferty et al. 2020
# Author: Paul Q. Sims
# Contact: paul.q.sims@gmail.com
# Date: 2020
################################################################################

# Readme - shows that the repository has been updated
rmarkdown::render("readme.Rmd", 
                  output_format = c("html_document", "github_document"))

# Data readme 
# rmarkdown::render("data/data_readme.Rmd", output_dir = "data/",
#                   output_format = c("html_document", "github_document"))

# Render .Rmd files of interest

docs_to_render <-
  c("R/analysis_summary.Rmd") # "R/analysis_repeatability.Rmd", "R/analysis_predictbehav.Rmd", "R/figures.Rmd", 

for (i in seq_along(docs_to_render)) {
  rmarkdown::render(docs_to_render[i], output_dir = "reports/",
                    output_format = "all")
}