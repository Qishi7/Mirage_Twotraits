local({
  local_lib <- file.path(getwd(), "Rlib")
  if (dir.exists(local_lib)) {
    .libPaths(c(local_lib, .libPaths()))
  }
})

## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}
