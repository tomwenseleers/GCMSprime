# Run from the repository root: Rscript GCMSprime/pkgdown/build.R
package <- if (file.exists("GCMSprime/DESCRIPTION")) "GCMSprime" else "."
assets <- file.path(package, "pkgdown/assets")
dir.create(assets, recursive = TRUE, showWarnings = FALSE)
logos <- c("logo.png", "logo-dark.svg", "logo-light.svg", "logo-dark.png", "logo-light.png")
stopifnot(all(file.copy(file.path(package, "man/figures", logos), assets, overwrite = TRUE)))
pkgdown::build_site(package, examples = FALSE, install = FALSE,
  new_process = FALSE, preview = FALSE)

destination <- pkgdown::as_pkgdown(package)$dst_path
# Remove the automatic home header after pkgdown inserts the package logo.
home <- file.path(destination, "index.html")
html <- xml2::read_html(home)
page_title <- yaml::read_yaml(file.path(package, "_pkgdown.yml"))$home$title
title_node <- xml2::xml_find_first(html, ".//head/title")
xml2::xml_text(title_node) <- page_title
xml2::xml_set_attr(xml2::xml_find_all(html,
  ".//head/meta[@property='og:title' or @name='twitter:title']"), "content", page_title)
xml2::xml_remove(xml2::xml_find_all(html, ".//main//div[@class='page-header']"), free = FALSE)
xml2::write_html(html, home)
pkgdown::build_search(package)
file.copy(file.path(package, "pkgdown/extra.css"),
  file.path(destination, "extra.css"), overwrite = TRUE)
file.copy(file.path(package, "pkgdown/favicon/site.webmanifest"),
  file.path(destination, "site.webmanifest"), overwrite = TRUE)
