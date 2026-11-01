package <- if(file.exists("GCMSprime/DESCRIPTION")) "GCMSprime" else "."
pkgdown::build_site(package,examples=FALSE,install=FALSE,new_process=FALSE,preview=FALSE)
destination<-pkgdown::as_pkgdown(package)$dst_path
home<-file.path(destination,"index.html")
html<-xml2::read_html(home)
xml2::xml_remove(xml2::xml_find_all(html,".//main//div[@class='page-header']"),free=FALSE)
xml2::write_html(html,home)
