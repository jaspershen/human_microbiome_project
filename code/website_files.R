# distill::create_website(dir = "website_files/", gh_pages = TRUE)
setwd(masstools::get_project_wd())
##render websites
rmarkdown::render_site(input = "website_files/")
##copy website docs to root folder
file.copy(from = "website_files/docs/", to = ".", recursive = TRUE, overwrite = TRUE)
