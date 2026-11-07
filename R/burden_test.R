gene_burden_test <- function(data, n_case, n_control, gene_col = "Gene",
                             case_col = "No.case", control_col = "No.contr",
                             variants_per_gene = NULL) {
  required <- c(gene_col, case_col, control_col)
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("data are missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  d <- data.frame(
    Gene = as.character(data[[gene_col]]),
    case = as.numeric(data[[case_col]]),
    control = as.numeric(data[[control_col]]),
    stringsAsFactors = FALSE
  )
  agg <- aggregate(cbind(case, control) ~ Gene, d, sum)
  if (is.null(variants_per_gene)) {
    variant_counts <- aggregate(case ~ Gene, d, length)
    names(variant_counts)[2] <- "variant_count"
    agg <- merge(agg, variant_counts, by = "Gene", all.x = TRUE)
  } else {
    agg$variant_count <- variants_per_gene
  }
  p_value <- odds_ratio <- numeric(nrow(agg))
  for (i in seq_len(nrow(agg))) {
    case_an <- 2 * n_case * agg$variant_count[i]
    control_an <- 2 * n_control * agg$variant_count[i]
    tab <- matrix(
      c(
        agg$case[i], case_an - agg$case[i],
        agg$control[i], control_an - agg$control[i]
      ),
      nrow = 2,
      byrow = TRUE
    )
    ft <- fisher.test(tab, alternative = "two.sided")
    p_value[i] <- ft$p.value
    odds_ratio[i] <- unname(ft$estimate)
  }
  agg$odds_ratio <- odds_ratio
  agg$p_value <- p_value
  agg$bh_fdr <- p.adjust(p_value, method = "BH")
  agg$score <- -log10(pmax(p_value, .Machine$double.xmin))
  agg
}
