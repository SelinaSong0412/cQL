#' Fit clustered Q-learning with the M-out-of-N cluster bootstrap
#'
#' Fits a two-stage clustered Q-learning analysis for one cSMART dataset using
#' the M-out-of-N cluster bootstrap described in the manuscript. The function
#' supports both:
#'
#' - type I cSMARTs, where every cluster is re-randomized at stage 2; and
#' - type II/III cSMARTs, where only a subset of clusters is re-randomized at
#'   stage 2 and the stage-2 treatment column is `NA` for clusters not
#'   re-randomized.
#'
#' The data must be at the individual level, with one row per individual.
#' Cluster-level treatments and stage-2 tailoring variables should be repeated
#' across individuals within each cluster.
#'
#' @param data A data frame containing one row per individual.
#' @param stage2_formula A formula for the stage-2 Q-function, for example
#'   `Y ~ X1 * A1 + A1 * A2 + A2:X2`.
#' @param stage1_formula A formula for the stage-1 Q-function, for example
#'   `Y1 ~ X1 * A1`. The response name does not need to exist in `data`; it is
#'   created internally as the stage-1 pseudo-outcome.
#' @param cluster Character scalar giving the cluster identifier column name.
#' @param stage1_treat Character scalar giving the stage-1 treatment column
#'   name.
#' @param stage2_treat Character scalar giving the stage-2 treatment column
#'   name. For limited second-stage re-randomization, clusters not
#'   re-randomized at stage 2 should have `NA` in this column.
#' @param stage2_tailoring_vars Character vector of cluster-level candidate
#'   tailoring variables that interact with `stage2_treat` in
#'   `stage2_formula`. For the example `Y ~ X1 * A1 + A1 * A2 + A2:X2`, use
#'   `c("A1", "X2")`.
#' @param alpha Significance level used for the bootstrap confidence intervals.
#'   Defaults to `0.05` for 95% confidence intervals.
#' @param n_boot Number of bootstrap replicates for each stage. Defaults to
#'   `1000`.
#' @param working_correlation Working correlation used when estimating the
#'   non-regularity measure that determines the stage-1 resample size.
#'   `"exchangeable"` is the manuscript default.
#' @param m_method Method for selecting the stage-1 cluster resample size.
#'   `"fixed_xi"` uses the manuscript's fixed tuning parameter rule,
#'   `"full_n"` uses the full-cluster bootstrap with `M = N`, and `"user"`
#'   uses the value supplied through `m`.
#' @param fixed_xi Tuning parameter used when `m_method = "fixed_xi"`. The
#'   manuscript default is `0.025`.
#' @param m User-supplied stage-1 cluster resample size. Required when
#'   `m_method = "user"`.
#' @param seed Optional random seed.
#' @param verbose Logical; if `TRUE`, prints the selected stage-1 resample size
#'   and the bootstrap progress warnings when relevant.
#'
#' @return An object of class `"cQL_fit"` containing:
#' - `stage_summary`: a data frame reporting each stage's total number of
#'   clusters `N`, number of randomized clusters `N_rand`, and bootstrap
#'   resample size `M`;
#' - `stage2`: a data frame with the stage-2 coefficient estimates,
#'   bootstrap standard errors, confidence intervals, p-values, and
#'   significance stars;
#' - `stage1`: the analogous stage-1 output;
#' - `p_hat`: the estimated degree of non-regularity;
#' - `stage1_m`: the selected stage-1 cluster resample size;
#' - `stage2_cluster_count`: the number of clusters re-randomized at stage 2;
#' - `cluster_count`: the total number of clusters;
#' - `design`: `"full_second_stage_rerandomization"` or
#'   `"limited_second_stage_rerandomization"`;
#' - `treatment_mapping`: the internal `-1/+1` recoding used for the two
#'   treatment columns.
#'
#' @examples
#' toy_data <- simulate_csmart_data(
#'   n_clusters = 48,
#'   cluster_size = 12,
#'   rerandomization = "nonresponder",
#'   p_rerand = 0.7,
#'   seed = 7
#' )
#'
#' fit <- cQL(
#'   data = toy_data,
#'   stage2_formula = Y ~ X1 * A1 + A1 * A2 + A2:X2,
#'   stage1_formula = Y1 ~ X1 * A1,
#'   cluster = "cluster_id",
#'   stage1_treat = "A1",
#'   stage2_treat = "A2",
#'   stage2_tailoring_vars = c("A1", "X2"),
#'   n_boot = 200,
#'   seed = 1,
#'   verbose = FALSE
#' )
#'
#' fit$stage_summary
#' fit
#' @export
cQL <- function(data,
                stage2_formula,
                stage1_formula,
                cluster,
                stage1_treat,
                stage2_treat,
                stage2_tailoring_vars,
                alpha = 0.05,
                n_boot = 1000,
                working_correlation = c("exchangeable", "independence"),
                m_method = c("fixed_xi", "full_n", "user"),
                fixed_xi = 0.025,
                m = NULL,
                seed = NULL,
                verbose = TRUE) {
  working_correlation <- match.arg(working_correlation)
  m_method <- match.arg(m_method)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  if (!is.character(cluster) || length(cluster) != 1L) {
    stop("`cluster` must be a single column name.", call. = FALSE)
  }
  if (!is.character(stage1_treat) || length(stage1_treat) != 1L) {
    stop("`stage1_treat` must be a single column name.", call. = FALSE)
  }
  if (!is.character(stage2_treat) || length(stage2_treat) != 1L) {
    stop("`stage2_treat` must be a single column name.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || n_boot < 50) {
    stop("`n_boot` must be a single number of at least 50.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  prepared <- .cql_prepare_data(
    data = data,
    stage2_formula = stage2_formula,
    stage1_formula = stage1_formula,
    cluster = cluster,
    stage1_treat = stage1_treat,
    stage2_treat = stage2_treat,
    stage2_tailoring_vars = stage2_tailoring_vars
  )

  stage2_fit <- stats::lm(stage2_formula, data = prepared$stage2_data)
  if (anyNA(stats::coef(stage2_fit))) {
    stop(
      "The stage-2 formula is not estimable on the provided data. ",
      "At least one stage-2 coefficient is NA.",
      call. = FALSE
    )
  }
  stage1_response <- .cql_formula_response(stage1_formula)
  outcome_name <- .cql_formula_response(stage2_formula)

  stage1_data <- prepared$data
  stage1_data[[stage1_response]] <- .cql_pseudo_outcome(
    stage2_fit = stage2_fit,
    data = prepared$data,
    stage2_treat = stage2_treat,
    stage2_available = prepared$data$.cql_stage2_available,
    outcome_name = outcome_name
  )
  stage1_fit <- stats::lm(stage1_formula, data = stage1_data)
  if (anyNA(stats::coef(stage1_fit))) {
    stop(
      "The stage-1 formula is not estimable on the provided data. ",
      "At least one stage-1 coefficient is NA.",
      call. = FALSE
    )
  }

  p_hat <- .cql_nonregularity(
    data = prepared$stage2_data,
    stage2_formula = stage2_formula,
    cluster = cluster,
    stage2_treat = stage2_treat,
    alpha = alpha,
    working_correlation = working_correlation
  )

  cluster_count <- length(unique(prepared$data[[cluster]]))
  stage2_cluster_count <- length(unique(prepared$stage2_data[[cluster]]))
  stage1_randomized_clusters <- cluster_count

  stage1_m <- .cql_choose_m(
    cluster_count = stage1_randomized_clusters,
    p_hat = p_hat,
    m_method = m_method,
    fixed_xi = fixed_xi,
    m = m
  )

  stage_summary <- .cql_stage_summary(
    total_clusters = cluster_count,
    stage1_randomized_clusters = stage1_randomized_clusters,
    stage2_randomized_clusters = stage2_cluster_count,
    stage1_m = stage1_m
  )

  if (isTRUE(verbose)) {
    message("Selected stage-1 resample size M = ", stage1_m)
  }

  stage2_star <- .cql_stage2_bootstrap(
    stage2_formula = stage2_formula,
    stage2_data = prepared$stage2_data,
    cluster = cluster,
    ref_coef = stats::coef(stage2_fit),
    n_boot = n_boot
  )

  stage1_star <- .cql_stage1_bootstrap(
    stage2_formula = stage2_formula,
    stage1_formula = stage1_formula,
    data = prepared$data,
    cluster = cluster,
    stage2_treat = stage2_treat,
    stage1_response = stage1_response,
    outcome_name = outcome_name,
    ref_coef_stage1 = stats::coef(stage1_fit),
    ref_coef_stage2 = stats::coef(stage2_fit),
    m = stage1_m,
    n_boot = n_boot
  )

  stage2_table <- .cql_bootstrap_table(
    estimate = stats::coef(stage2_fit),
    star_draws = stage2_star,
    alpha = alpha
  )
  stage1_table <- .cql_bootstrap_table(
    estimate = stats::coef(stage1_fit),
    star_draws = stage1_star,
    alpha = alpha
  )

  fit <- list(
    call = match.call(),
    stage_summary = stage_summary,
    stage2 = stage2_table,
    stage1 = stage1_table,
    p_hat = p_hat,
    stage1_m = stage1_m,
    cluster_count = cluster_count,
    stage2_cluster_count = stage2_cluster_count,
    design = if (all(prepared$data$.cql_stage2_available)) {
      "full_second_stage_rerandomization"
    } else {
      "limited_second_stage_rerandomization"
    },
    treatment_mapping = prepared$treatment_mapping,
    working_correlation = working_correlation,
    fixed_xi = fixed_xi,
    alpha = alpha
  )

  class(fit) <- "cQL_fit"
  fit
}

#' @export
print.cQL_fit <- function(x, ...) {
  cat("cQL fit\n")
  cat("Design:", x$design, "\n")
  cat("Estimated non-regularity p_hat:", round(x$p_hat, 4), "\n")
  cat("\nStage summary\n")
  print(x$stage_summary, row.names = FALSE)
  cat("\n")

  cat("Stage 2\n")
  print(x$stage2, row.names = FALSE)
  cat("\nStage 1\n")
  print(x$stage1, row.names = FALSE)
  invisible(x)
}

.cql_stage_summary <- function(total_clusters,
                               stage1_randomized_clusters,
                               stage2_randomized_clusters,
                               stage1_m) {
  data.frame(
    stage = c("Stage 2", "Stage 1"),
    N = c(total_clusters, total_clusters),
    N_rand = c(stage2_randomized_clusters, stage1_randomized_clusters),
    M = c(stage2_randomized_clusters, stage1_m),
    bootstrap = c(
      "Cluster bootstrap on stage-2 randomized clusters",
      "M-out-of-N cluster bootstrap"
    ),
    row.names = NULL,
    check.names = FALSE
  )
}

.cql_prepare_data <- function(data,
                              stage2_formula,
                              stage1_formula,
                              cluster,
                              stage1_treat,
                              stage2_treat,
                              stage2_tailoring_vars) {
  needed <- unique(c(
    cluster,
    stage1_treat,
    stage2_treat,
    all.vars(stage2_formula),
    all.vars(stage1_formula)
  ))
  needed <- needed[needed != .cql_formula_response(stage1_formula)]
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      "The following required columns are missing from `data`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  .cql_validate_stage2_formula(
    stage2_formula = stage2_formula,
    stage2_treat = stage2_treat,
    stage2_tailoring_vars = stage2_tailoring_vars
  )

  out <- data
  mapping1 <- .cql_recode_binary_treatment(out[[stage1_treat]], stage1_treat)
  out[[stage1_treat]] <- mapping1$values

  stage2_available <- !is.na(out[[stage2_treat]])
  if (!any(stage2_available)) {
    stop("No stage-2 randomized clusters were found in `data`.", call. = FALSE)
  }
  mapping2 <- .cql_recode_binary_treatment(
    out[[stage2_treat]],
    stage2_treat,
    allow_na = TRUE
  )
  out[[stage2_treat]] <- mapping2$values
  out$.cql_stage2_available <- stage2_available

  response_stage2 <- .cql_formula_response(stage2_formula)
  response_stage1 <- .cql_formula_response(stage1_formula)

  required_vars <- unique(c(
    cluster,
    stage1_treat,
    setdiff(all.vars(stage2_formula), stage2_treat),
    setdiff(all.vars(stage1_formula), response_stage1)
  ))
  keep_complete <- stats::complete.cases(out[, required_vars, drop = FALSE])
  out <- out[keep_complete, , drop = FALSE]
  out$.cql_stage2_available <- !is.na(out[[stage2_treat]])

  if (anyNA(out[[response_stage2]])) {
    stop("The stage-2 outcome cannot contain missing values.", call. = FALSE)
  }

  .cql_check_cluster_consistency(out, cluster, stage1_treat, allow_na = FALSE)
  .cql_check_cluster_consistency(out, cluster, stage2_treat, allow_na = TRUE)
  for (var in stage2_tailoring_vars) {
    .cql_check_cluster_consistency(out, cluster, var, allow_na = FALSE)
  }

  if (!response_stage1 %in% names(out)) {
    out[[response_stage1]] <- 0
  }

  stage2_data <- out[out$.cql_stage2_available, , drop = FALSE]
  if (length(unique(stage2_data[[cluster]])) < 2L) {
    stop(
      "At least two stage-2 randomized clusters are required.",
      call. = FALSE
    )
  }

  list(
    data = out,
    stage2_data = stage2_data,
    treatment_mapping = list(
      stage1 = mapping1$mapping,
      stage2 = mapping2$mapping
    )
  )
}

.cql_validate_stage2_formula <- function(stage2_formula,
                                         stage2_treat,
                                         stage2_tailoring_vars) {
  term_labels <- attr(stats::terms(stage2_formula), "term.labels")
  treat_regex <- paste0("(^|:)", stage2_treat, "(:|$)")
  treat_terms <- term_labels[grepl(treat_regex, term_labels)]

  if (!any(treat_terms == stage2_treat)) {
    stop(
      "The stage-2 formula must include the main effect for `stage2_treat`.",
      call. = FALSE
    )
  }

  interaction_terms <- setdiff(treat_terms, stage2_treat)
  for (term in interaction_terms) {
    pieces <- strsplit(term, ":", fixed = TRUE)[[1]]
    other <- setdiff(pieces, stage2_treat)
    if (length(other) != 1L || !(other %in% stage2_tailoring_vars)) {
      stop(
        "Every interaction involving `stage2_treat` must be a two-way ",
        "interaction with one of `stage2_tailoring_vars`.",
        call. = FALSE
      )
    }
  }
}

.cql_formula_response <- function(formula_obj) {
  vars <- all.vars(formula_obj)
  if (length(vars) == 0L) {
    stop("Formulas must include a response.", call. = FALSE)
  }
  vars[[1]]
}

.cql_recode_binary_treatment <- function(x, name, allow_na = FALSE) {
  observed <- x[!is.na(x)]
  levels_obs <- unique(observed)
  if (length(levels_obs) != 2L) {
    stop(
      "Column `", name, "` must contain exactly two non-missing treatment ",
      "levels.",
      call. = FALSE
    )
  }

  if (is.factor(x)) {
    level_order <- levels(x)[levels(x) %in% levels_obs]
  } else if (is.numeric(observed) && setequal(sort(unique(observed)), c(-1, 1))) {
    level_order <- c(-1, 1)
  } else {
    level_order <- sort(unique(as.character(observed)))
  }

  if (length(level_order) != 2L) {
    level_order <- as.character(levels_obs)
  }

  recoded <- rep(NA_real_, length(x))
  recoded[!is.na(x) & as.character(x) == as.character(level_order[1])] <- -1
  recoded[!is.na(x) & as.character(x) == as.character(level_order[2])] <- 1

  if (!allow_na && anyNA(recoded)) {
    stop("Column `", name, "` contains missing treatment values.", call. = FALSE)
  }

  list(
    values = recoded,
    mapping = data.frame(
      original = as.character(level_order),
      coded = c(-1, 1),
      row.names = NULL
    )
  )
}

.cql_check_cluster_consistency <- function(data, cluster, var, allow_na) {
  split_vals <- split(data[[var]], data[[cluster]])
  bad <- vapply(
    split_vals,
    function(x) {
      x <- unique(x)
      x <- x[!is.na(x)]
      if (allow_na && length(x) == 0L) {
        return(FALSE)
      }
      length(x) > 1L
    },
    logical(1)
  )
  if (any(bad)) {
    stop(
      "Column `", var, "` must be constant within each cluster.",
      call. = FALSE
    )
  }
}

.cql_pseudo_outcome <- function(stage2_fit,
                                data,
                                stage2_treat,
                                stage2_available,
                                outcome_name) {
  plus_data <- data
  minus_data <- data
  plus_data[[stage2_treat]] <- 1
  minus_data[[stage2_treat]] <- -1

  q_plus <- stats::predict(stage2_fit, newdata = plus_data)
  q_minus <- stats::predict(stage2_fit, newdata = minus_data)
  q_opt <- pmax(q_plus, q_minus)

  if (all(stage2_available)) {
    return(q_opt)
  }

  pseudo <- data[[outcome_name]]
  pseudo[stage2_available] <- q_opt[stage2_available]
  pseudo
}

.cql_nonregularity <- function(data,
                               stage2_formula,
                               cluster,
                               stage2_treat,
                               alpha,
                               working_correlation) {
  stage2_fit_lm <- stats::lm(stage2_formula, data = data)
  working_fit <- try(
    geepack::geeglm(
      formula = stage2_formula,
      data = data,
      id = data[[cluster]],
      corstr = working_correlation
    ),
    silent = TRUE
  )

  if (inherits(working_fit, "try-error")) {
    coef_work <- stats::coef(stage2_fit_lm)
    cov_work <- sandwich::vcovCL(stage2_fit_lm, cluster = data[[cluster]])
  } else {
    coef_work <- stats::coef(working_fit)
    cov_work <- summary(working_fit)$cov.scaled
  }

  keep_coef <- !is.na(coef_work)
  coef_work <- coef_work[keep_coef]
  coef_names <- names(coef_work)
  cov_names <- intersect(coef_names, rownames(cov_work))
  coef_work <- coef_work[cov_names]
  cov_work <- cov_work[cov_names, cov_names, drop = FALSE]

  cluster_rows <- !duplicated(data[[cluster]])
  cluster_data <- data[cluster_rows, , drop = FALSE]
  plus_data <- cluster_data
  minus_data <- cluster_data
  plus_data[[stage2_treat]] <- 1
  minus_data[[stage2_treat]] <- -1

  x_plus <- stats::model.matrix(stage2_formula, data = plus_data)
  x_minus <- stats::model.matrix(stage2_formula, data = minus_data)
  contrast_design <- (x_plus - x_minus) / 2
  contrast_design <- contrast_design[, names(coef_work), drop = FALSE]

  contrast_est <- drop(contrast_design %*% coef_work)
  contrast_var <- rowSums((contrast_design %*% cov_work) * contrast_design)
  contrast_se <- sqrt(pmax(contrast_var, 0))

  t_stat <- mapply(
    .cql_safe_t_stat,
    estimate = contrast_est,
    se = contrast_se,
    USE.NAMES = FALSE
  )

  cluster_sizes <- as.numeric(table(data[[cluster]]))
  eta <- stats::qt(
    p = 1 - alpha / (2 * nrow(cluster_data)),
    df = pmax(cluster_sizes - 1, 1)
  )

  mean(abs(t_stat) <= eta)
}

.cql_safe_t_stat <- function(estimate, se) {
  if (is.na(se) || se <= .Machine$double.eps) {
    if (abs(estimate) <= .Machine$double.eps) {
      return(0)
    }
    return(Inf)
  }
  estimate / se
}

.cql_choose_m <- function(cluster_count, p_hat, m_method, fixed_xi, m) {
  if (m_method == "full_n") {
    return(cluster_count)
  }

  if (m_method == "user") {
    if (is.null(m) || length(m) != 1L || m < 1 || m > cluster_count) {
      stop(
        "When `m_method = \"user\"`, `m` must be an integer between 1 and N.",
        call. = FALSE
      )
    }
    return(as.integer(ceiling(m)))
  }

  if (!is.numeric(fixed_xi) || length(fixed_xi) != 1L || fixed_xi <= 0) {
    stop("`fixed_xi` must be a positive number.", call. = FALSE)
  }

  exponent <- 1 - p_hat * (fixed_xi / (1 + fixed_xi))
  as.integer(min(cluster_count, max(1L, ceiling(cluster_count^exponent))))
}

.cql_stage2_bootstrap <- function(stage2_formula,
                                  stage2_data,
                                  cluster,
                                  ref_coef,
                                  n_boot) {
  .cql_bootstrap_loop(
    n_boot = n_boot,
    ref_coef = ref_coef,
    max_iter = 10L * n_boot,
    worker = function() {
      boot_data <- .cql_cluster_bootstrap_sample(
        data = stage2_data,
        cluster = cluster,
        n_clusters = length(unique(stage2_data[[cluster]]))
      )
      fit <- try(stats::lm(stage2_formula, data = boot_data), silent = TRUE)
      if (inherits(fit, "try-error")) {
        return(NULL)
      }
      coef_boot <- .cql_align_coef(ref_coef, stats::coef(fit))
      if (anyNA(coef_boot)) {
        return(NULL)
      }
      2 * ref_coef - coef_boot
    }
  )
}

.cql_stage1_bootstrap <- function(stage2_formula,
                                  stage1_formula,
                                  data,
                                  cluster,
                                  stage2_treat,
                                  stage1_response,
                                  outcome_name,
                                  ref_coef_stage1,
                                  ref_coef_stage2,
                                  m,
                                  n_boot) {
  .cql_bootstrap_loop(
    n_boot = n_boot,
    ref_coef = ref_coef_stage1,
    max_iter = 15L * n_boot,
    worker = function() {
      boot_data <- .cql_cluster_bootstrap_sample(
        data = data,
        cluster = cluster,
        n_clusters = m
      )
      boot_stage2 <- boot_data[boot_data$.cql_stage2_available, , drop = FALSE]
      if (length(unique(boot_stage2[[cluster]])) < 2L) {
        return(NULL)
      }

      fit_stage2 <- try(stats::lm(stage2_formula, data = boot_stage2), silent = TRUE)
      if (inherits(fit_stage2, "try-error")) {
        return(NULL)
      }
      coef_stage2 <- .cql_align_coef(ref_coef_stage2, stats::coef(fit_stage2))
      if (anyNA(coef_stage2)) {
        return(NULL)
      }

      boot_data[[stage1_response]] <- .cql_pseudo_outcome(
        stage2_fit = fit_stage2,
        data = boot_data,
        stage2_treat = stage2_treat,
        stage2_available = boot_data$.cql_stage2_available,
        outcome_name = outcome_name
      )

      fit_stage1 <- try(stats::lm(stage1_formula, data = boot_data), silent = TRUE)
      if (inherits(fit_stage1, "try-error")) {
        return(NULL)
      }
      coef_stage1 <- .cql_align_coef(ref_coef_stage1, stats::coef(fit_stage1))
      if (anyNA(coef_stage1)) {
        return(NULL)
      }

      2 * ref_coef_stage1 - coef_stage1
    }
  )
}

.cql_bootstrap_loop <- function(n_boot, ref_coef, max_iter, worker) {
  out <- matrix(
    NA_real_,
    nrow = length(ref_coef),
    ncol = n_boot,
    dimnames = list(names(ref_coef), NULL)
  )

  success <- 0L
  iter <- 0L
  while (success < n_boot && iter < max_iter) {
    iter <- iter + 1L
    draw <- worker()
    if (is.null(draw)) {
      next
    }
    success <- success + 1L
    out[, success] <- draw
  }

  if (success == 0L) {
    stop("All bootstrap refits failed.", call. = FALSE)
  }
  if (success < n_boot) {
    warning(
      "Only ", success, " successful bootstrap refits were obtained out of ",
      n_boot, ".",
      call. = FALSE
    )
    out <- out[, seq_len(success), drop = FALSE]
  }

  out
}

.cql_cluster_bootstrap_sample <- function(data, cluster, n_clusters) {
  cluster_ids <- unique(data[[cluster]])
  sampled_ids <- sample(cluster_ids, size = n_clusters, replace = TRUE)
  rows <- unlist(
    lapply(sampled_ids, function(id) which(data[[cluster]] == id)),
    use.names = FALSE
  )
  data[rows, , drop = FALSE]
}

.cql_align_coef <- function(reference, candidate) {
  out <- rep(NA_real_, length(reference))
  names(out) <- names(reference)
  out[names(candidate)] <- candidate
  out
}

.cql_bootstrap_table <- function(estimate, star_draws, alpha) {
  lower <- apply(
    star_draws,
    1,
    stats::quantile,
    probs = alpha / 2,
    na.rm = TRUE,
    names = FALSE
  )
  upper <- apply(
    star_draws,
    1,
    stats::quantile,
    probs = 1 - alpha / 2,
    na.rm = TRUE,
    names = FALSE
  )
  se <- apply(star_draws, 1, stats::sd, na.rm = TRUE)
  p_value <- apply(star_draws, 1, .cql_boot_p_value)

  data.frame(
    term = names(estimate),
    estimate = unname(estimate),
    conf.low = unname(lower),
    conf.high = unname(upper),
    std.error = unname(se),
    p.value = unname(p_value),
    significance = .cql_significance_stars(p_value),
    row.names = NULL,
    check.names = FALSE
  )
}

.cql_boot_p_value <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0L) {
    return(NA_real_)
  }
  left <- (sum(x <= 0) + 1) / (n + 1)
  right <- (sum(x >= 0) + 1) / (n + 1)
  min(1, 2 * min(left, right))
}

.cql_significance_stars <- function(p_value) {
  stars <- rep("", length(p_value))
  stars[p_value < 0.05] <- "*"
  stars[p_value < 0.01] <- "**"
  stars[p_value < 0.001] <- "***"
  stars
}
