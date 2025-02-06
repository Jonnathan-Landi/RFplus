RFplusNews <- function() {
  newsfile <- file.path(system.file(package="RFplus"), "NEWS")
  file.show(newsfile)
}

