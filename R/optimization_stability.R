# Optional maximum-likelihood refinement. The EM implementation and all
# variant Bayes factors are deliberately independent of this file.

# Concave maximization in the ORIGINAL four simplex weights. Newton steps
# remain on the current face; pairwise vertex exchanges permit zero weights
# to re-enter. No positive lower bound or statistical penalty is introduced.
mirage_simplex_optimize <- function(log_components, pi, tolerance = 1e-7,
                                    max_iter = 2000L) {
  L <- as.matrix(log_components)
  if (ncol(L) != 4L || !nrow(L) || any(!is.finite(L)) ||
      length(pi) != 4L || any(!is.finite(pi)) || any(pi < 0) ||
      abs(sum(pi) - 1) > 1e-10 || !is.finite(tolerance) || tolerance <= 0 ||
      max_iter < 1) stop("Invalid simplex optimization inputs.", call. = FALSE)
  p <- as.numeric(pi) / sum(pi)
  # Scaling each gene's components leaves the objective and derivatives intact.
  shift <- apply(L, 1L, max)
  F <- exp(L - shift)
  evaluate <- function(p) {
    den <- drop(F %*% p)
    if (any(den <= 0)) stop("Component dynamic range exceeds numeric precision.", call. = FALSE)
    R <- F / den
    list(ll = sum(log(den) + shift), R = R, g = colSums(R))
  }
  z <- evaluate(p)
  history <- z$ll
  reason <- "iteration_limit"
  for (it in seq_len(max_iter)) {
    gap <- max(0, max(z$g) - sum(p * z$g))
    if (gap <= tolerance) { reason <- "simplex_gap"; break }
    active <- which(p > 0)
    direction <- numeric(4)
    # Reduced Newton system on the active simplex face, without a chart.
    if (length(active) > 1L) {
      base <- tail(active, 1L); free <- head(active, -1L)
      X <- z$R[, free, drop = FALSE] - z$R[, base]
      h <- crossprod(X)
      grad <- z$g[free] - z$g[base]
      step <- tryCatch(drop(solve(h, grad)), error = function(e) NULL)
      if (!is.null(step) && all(is.finite(step))) {
        direction[free] <- step; direction[base] <- -sum(step)
      }
    }
    slope <- sum(z$g * direction)
    # Use a pairwise Frank-Wolfe/away step if the active face is stationary
    # or the Newton direction has insufficient ascent. This escapes vertices.
    if (!is.finite(slope) || slope <= 0 ||
        diff(range(z$g[active])) < gap * .1) {
      direction[] <- 0
      direction[which.max(z$g)] <- 1
      away <- active[which.min(z$g[active])]
      direction[away] <- direction[away] - 1
    }
    neg <- which(direction < 0)
    limit <- min(p[neg] / -direction[neg])
    rd <- drop(z$R %*% direction)
    derivative <- function(a) sum(rd / (1 + a * rd))
    # Exact one-dimensional concave maximization, including a boundary step.
    edge <- derivative(limit)
    alpha <- if (!is.na(edge) && edge >= 0) limit else
      stats::uniroot(derivative, c(0, limit), tol = .Machine$double.eps^0.75)$root
    pn <- p + alpha * direction
    if (alpha == limit) pn[neg[which.min(p[neg] / -direction[neg])]] <- 0
    pn <- pmax(pn, 0); pn <- pn / sum(pn)
    zn <- evaluate(pn)
    if (zn$ll < z$ll - 1e-10 * (1 + abs(z$ll))) {
      reason <- "numerical_decrease"; break
    }
    if (max(abs(pn - p)) == 0) { reason <- "numeric_stall"; break }
    p <- pn; z <- zn; history <- c(history, z$ll)
  }
  gap <- max(0, max(z$g) - sum(p * z$g))
  names(p) <- c("00", "10", "01", "11")
  list(pi = p, loglik = z$ll, gradient = z$g, fw_gap = gap,
       # By concavity, gap bounds the remaining fixed-eta objective error.
       certified = gap <= tolerance, reason = reason,
       iterations = length(history) - 1L, loglik_trace = history)
}

# Compress identical nonzero variant BFs only; multiplicities are exact.
# Zero BF rows have exactly zero log-mixture and eta derivative.
.mirage_optimization_trait <- function(d, n, categories) {
  d <- d[d$log_bf != 0, , drop = FALSE]
  if (!nrow(d)) return(list(gene = integer(), category = integer(),
                           lb = numeric(), count = numeric()))
  cls <- match(d$category, categories)
  key <- paste(d$gene_id, cls, sprintf("%.17g", d$log_bf), sep = "|")
  keep <- !duplicated(key)
  ix <- match(key, key[keep])
  list(gene = d$gene_id[keep], category = cls[keep], lb = d$log_bf[keep],
       count = tabulate(ix, nbins = sum(keep)))
}

.mirage_optimization_model <- function(d1, d2, n, k1, k2) {
  trait <- function(d, eta, k) {
    L <- numeric(n); U <- matrix(0, n, k)
    if (length(d$gene)) {
      e <- eta[d$category]
      # Stable mixture and derivative, including exact eta endpoints.
      a <- log1p(-e); b <- log(e) + d$lb; m <- pmax(a, b)
      lm <- m + log(exp(a - m) + exp(b - m))
      der <- numeric(length(e)); pos <- d$lb >= 0
      der[pos] <- exp(d$lb[pos] - lm[pos]) * (-expm1(-d$lb[pos]))
      der[!pos] <- exp(-lm[!pos]) * expm1(d$lb[!pos])
      sums <- rowsum(cbind(d$count * lm), d$gene, reorder = FALSE)
      L[as.integer(rownames(sums))] <- sums[, 1L]
      sums <- rowsum(cbind(d$count * der), d$gene + n * (d$category - 1L), reorder = FALSE)
      U[as.integer(rownames(sums))] <- sums[, 1L]
    }
    list(L = L, U = U)
  }
  components <- function(eta) {
    a <- trait(d1, eta[seq_len(k1)], k1)
    b <- trait(d2, eta[k1 + seq_len(k2)], k2)
    list(L = cbind(0, a$L, b$L, a$L + b$L), u = a$U, v = b$U)
  }
  evaluate <- function(pi, eta, components_cached = NULL) {
    c <- if (is.null(components_cached)) components(eta) else components_cached
    W <- sweep(c$L, 2L, log(pi), "+")
    m <- apply(W, 1L, max)
    den <- m + log(rowSums(exp(W - m)))
    tau <- exp(W - den)
    list(ll = sum(den), grad_pi = colSums(exp(c$L - den)),
         grad_eta = c(drop(crossprod(c$u, tau[, 2L] + tau[, 4L])),
                      drop(crossprod(c$v, tau[, 3L] + tau[, 4L]))),
         tau = tau, components = c$L)
  }
  list(n = n, k1 = k1, k2 = k2, components = components, evaluate = evaluate)
}

.mirage_refine_model <- function(model, em, estimate_pi, estimate_eta,
                                 starts = list(), tolerance_pi = 1e-7,
                                 tolerance_eta = 1e-6, max_iter = 3000L) {
  validate <- function(s) {
    p <- as.numeric(s$pi); e <- as.numeric(s$eta)
    if (length(p) != 4L || length(e) != model$k1 + model$k2 ||
        any(!is.finite(c(p, e))) || any(p < 0) || abs(sum(p) - 1) > 1e-10 ||
        any(e < 0 | e > 1)) stop("Invalid refinement start.", call. = FALSE)
    if (!estimate_pi && max(abs(p - em$pi)) > 1e-12)
      stop("Starts must preserve fixed pi.", call. = FALSE)
    if (!estimate_eta && max(abs(e - em$eta)) > 1e-12)
      stop("Starts must preserve fixed eta.", call. = FALSE)
    list(pi = p / sum(p), eta = e)
  }
  em <- validate(em)
  if (estimate_eta) {
    # Generic deterministic starts, independent of simulation truth and AUC.
    defaults <- list(
      list(pi = if (estimate_pi) rep(.25, 4) else em$pi, eta = rep(.1, length(em$eta))),
      list(pi = if (estimate_pi) c(.85, .05, .05, .05) else em$pi,
           eta = rep(.5, length(em$eta))))
  } else defaults <- list()
  candidates <- c(list(EM = em), defaults, starts)
  names(candidates) <- paste0("start", seq_along(candidates))
  names(candidates)[1L] <- "EM"
  records <- list()
  # Retain the untouched EM solution, even if every attempted optimizer fails.
  retain <- function(s, label, code, message) {
    ev <- model$evaluate(s$pi, s$eta)
    ge <- ev$grad_eta
    ge[s$eta == 0 & ge < 0 | s$eta == 1 & ge > 0] <- 0
    gap <- max(0, max(ev$grad_pi) - sum(s$pi * ev$grad_pi))
    list(pi = s$pi, eta = s$eta, ll = ev$ll, start = label,
         code = code, message = message, fw_gap = gap,
         eta_kkt = max(abs(ge)),
         stationary = (!estimate_pi || gap <= tolerance_pi) &&
                      (!estimate_eta || max(abs(ge)) <= tolerance_eta))
  }
  records[[1L]] <- retain(em, "EM_retained", NA_integer_, "Original EM retained unchanged")
  for (j in seq_along(candidates)) {
    s <- validate(candidates[[j]])
    initial <- retain(s, paste0(names(candidates)[j], "_initial"), NA_integer_,
                      "Feasible start retained unchanged")
    records[[length(records) + 1L]] <- initial
    best_seen <- initial
    attempt <- tryCatch({
      # Profile out concave simplex weights. The envelope gradient in eta is
      # evaluated at the certified pi solution; no stick-breaking variables.
      last <- NULL
      profile <- function(eta) {
        if (!is.null(last) && identical(eta, last$eta)) return(last)
        cc <- model$components(eta)
        po <- if (estimate_pi) mirage_simplex_optimize(cc$L, s$pi,
                    tolerance = tolerance_pi) else list(pi = s$pi, certified = TRUE)
        seen <- retain(list(pi = po$pi, eta = eta),
                       paste0(names(candidates)[j], "_evaluated"), NA_integer_,
                       "Best feasible objective evaluated along this start")
        if (seen$ll > best_seen$ll) best_seen <<- seen
        if (!po$certified) stop("Pi subproblem did not meet its FW certificate.")
        z <- model$evaluate(po$pi, eta, cc)
        last <<- c(list(pi = po$pi, eta = eta), z[c("ll", "grad_eta")])
        last
      }
      if (estimate_eta) {
        op <- stats::optim(s$eta, function(e) -profile(e)$ll,
                          function(e) -profile(e)$grad_eta, method = "L-BFGS-B",
                          lower = 0, upper = 1,
                          control = list(maxit = max_iter, factr = 1,
                                         pgtol = tolerance_eta))
        z <- profile(op$par)
        retain(z, names(candidates)[j], op$convergence,
               if (is.null(op$message)) "" else op$message)
      } else {
        z <- profile(s$eta)
        retain(z, names(candidates)[j], 0L, "Fixed eta; simplex optimization")
      }
    }, error = function(e) list(error = conditionMessage(e), start = names(candidates)[j]))
    records[[length(records) + 1L]] <- attempt
    records[[length(records) + 1L]] <- best_seen
  }
  valid <- which(vapply(records, function(z) !is.null(z$ll) && is.finite(z$ll), logical(1)))
  best <- records[[valid[which.max(vapply(records[valid], `[[`, numeric(1), "ll"))]]]
  best$estimate_pi <- estimate_pi; best$estimate_eta <- estimate_eta
  best$tolerance_pi <- tolerance_pi; best$tolerance_eta <- tolerance_eta
  best$em_loglik <- records[[1L]]$ll
  best$starts <- records
  best
}

# Explicit flags are required: callers must preserve their analysis's fixed
# parameters. This returns a separate result and never mutates the EM fit.
mirage_refine <- function(fit, estimate_pi, estimate_eta, starts = list(),
                          tolerance_pi = 1e-7, tolerance_eta = 1e-6,
                          max_iter = 3000L) {
  if (!inherits(fit, "mirage_two_trait")) stop("Expected a mirage_two_trait EM fit.")
  if (!is.null(fit$call$dirichlet_alpha))
    stop("Refinement supports unpenalized ML only; MAP fits are not supported.")
  if (length(estimate_pi) != 1 || !is.logical(estimate_pi) || is.na(estimate_pi) ||
      length(estimate_eta) != 1 || !is.logical(estimate_eta) || is.na(estimate_eta))
    stop("Supply explicit logical estimate_pi and estimate_eta flags.")
  genes <- fit$posterior$Gene
  c1 <- names(fit$parameters$eta_trait1); c2 <- names(fit$parameters$eta_trait2)
  model <- .mirage_optimization_model(
    .mirage_optimization_trait(fit$variant_info$trait1, length(genes), c1),
    .mirage_optimization_trait(fit$variant_info$trait2, length(genes), c2),
    length(genes), length(c1), length(c2))
  em <- list(pi = fit$parameters$pi,
             eta = c(fit$parameters$eta_trait1, fit$parameters$eta_trait2))
  ans <- .mirage_refine_model(model, em, estimate_pi, estimate_eta, starts,
                             tolerance_pi, tolerance_eta, max_iter)
  z <- model$evaluate(ans$pi, ans$eta)
  ans$parameters <- list(pi = stats::setNames(ans$pi, c("00", "10", "01", "11")),
    eta_trait1 = stats::setNames(ans$eta[seq_along(c1)], c1),
    eta_trait2 = stats::setNames(ans$eta[length(c1) + seq_along(c2)], c2))
  ans$posterior <- data.frame(Gene = genes, z$tau,
    PP_trait1 = z$tau[, 2] + z$tau[, 4], PP_trait2 = z$tau[, 3] + z$tau[, 4],
    PP_pleiotropy = z$tau[, 4], check.names = FALSE)
  names(ans$posterior)[2:5] <- paste0("tau_", c("00", "10", "01", "11"))
  ans$bfdr <- list(trait1 = bfdr_table(genes, ans$posterior$PP_trait1, fit$trait_names[1]),
                   trait2 = bfdr_table(genes, ans$posterior$PP_trait2, fit$trait_names[2]),
                   pleiotropy = bfdr_table(genes, ans$posterior$PP_pleiotropy, "pleiotropy"))
  class(ans) <- "mirage_optimization"
  ans
}
