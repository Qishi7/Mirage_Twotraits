binary_auc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || all(truth)) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  n1 <- sum(truth)
  n0 <- sum(!truth)
  (sum(ranks[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
