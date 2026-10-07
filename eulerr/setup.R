#!/usr/bin/env Rscript
# Install exactly eulerr 8.3.1 in this folder. Run once; reproduction is offline.
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
local_library <- file.path(root, ".library")
dir.create(local_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(local_library, .libPaths()))
required_version <- "8.3.1"
stopifnot(getRversion() >= "4.2.0")

installed <- requireNamespace("eulerr", quietly = TRUE) &&
  as.character(packageVersion("eulerr")) == required_version
if (!installed) {
  repos <- "https://cloud.r-project.org"
  options(timeout = max(300, getOption("timeout")))
  # Use an exact-version binary when available; otherwise use its source archive.
  types <- unique(c(.Platform$pkgType, "source"))
  selected_url <- NULL
  selected_type <- "source"
  for (type in types) {
    available <- tryCatch(
      available.packages(repos = repos, type = type),
      error = function(e) NULL
    )
    if (!is.null(available) && "eulerr" %in% rownames(available) &&
        available["eulerr", "Version"] == required_version) {
      extension <- if (type == "source") ".tar.gz" else if (type == "win.binary") ".zip" else ".tgz"
      selected_url <- paste0(contrib.url(repos, type), "/eulerr_", required_version, extension)
      selected_type <- type
      break
    }
  }
  if (is.null(selected_url)) {
    selected_url <- paste0(repos, "/src/contrib/Archive/eulerr/eulerr_", required_version, ".tar.gz")
  }
  if (selected_type == "source") {
    message("A source build requires Rust >= 1.88.0 and the usual R build tools.")
  }
  install.packages(selected_url, repos = NULL, type = selected_type, lib = local_library)
}
stopifnot(as.character(packageVersion("eulerr", lib.loc = .libPaths())) == required_version)
message("Ready: eulerr ", required_version)
