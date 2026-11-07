core_test_data <- function(n_genes = 80, variants_per_gene = 6, ...) {
  n <- n_genes * variants_per_gene
  make <- function(phase) data.frame(
    Gene = rep(sprintf("G%03d", seq_len(n_genes)), each = variants_per_gene),
    ID = paste0("v", seq_len(n)), No.case = 0, No.contr = 0,
    category = rep(c("1", "2"), length.out = n),
    log_bf = 2 * sin(seq_len(n) / 5 + phase) - 0.5)
  list(trait1 = make(0), trait2 = make(.7))
}

toy_fit <- function(...) {
  d <- core_test_data()
  list(fit = mirage_two_trait(d$trait1, d$trait2, 500, 500, 500, 500,
    log_bf_col_trait1 = "log_bf", log_bf_col_trait2 = "log_bf",
    eta_init_trait1 = c("1" = .2, "2" = .05),
    eta_init_trait2 = c("1" = .2, "2" = .05), max_iter = 100, ...))
}
