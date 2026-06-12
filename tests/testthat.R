Sys.setenv(R_LIBS_USER = file.path(getwd(), "Rlib"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
testthat::test_dir("tests/testthat", reporter = "summary")

