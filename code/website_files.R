# distill::create_website(dir = "website_files/", gh_pages = TRUE)
<<<<<<< HEAD
masstools::setwd_project()
##render websites
rmarkdown::render_site(input = "website_files/")
##copy website docs to root folder
=======

tinyTools::setwd_project()
rmarkdown::render_site(input = "website_files/")
>>>>>>> 323e9fbbc4000cd0172a396a80c47cde0475b4da
file.copy(from = "website_files/docs/", to = ".", recursive = TRUE, overwrite = TRUE)
