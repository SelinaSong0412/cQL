#' Simulate a toy two-stage cSMART dataset
#'
#' Generates a simple individual-level cSMART dataset for examples and package
#' documentation. The returned data are in the format expected by [cQL()]:
#' one row per individual, a cluster identifier, separate stage-1 and stage-2
#' treatment columns, and `NA` in the stage-2 treatment column for clusters not
#' re-randomized at stage 2.
#'
#' @param n_clusters Number of clusters.
#' @param cluster_size Either a single common cluster size or an integer vector
#'   of length `n_clusters`.
#' @param rerandomization Either `"full"` for type I cSMART-style data, where
#'   all clusters are re-randomized at stage 2, `"nonresponder"` for a
#'   limited-re-randomization setting in which non-responder clusters are
#'   re-randomized, or `"responder"` for a limited-re-randomization setting in
#'   which responder clusters are re-randomized.
#' @param p_rerand Approximate proportion of clusters re-randomized at stage 2
#'   when `rerandomization` is `"nonresponder"` or `"responder"`.
#' @param icc Approximate intra-cluster correlation for the individual-level
#'   outcome.
#' @param seed Optional random seed.
#'
#' @return A data frame with one row per individual.
#'
#' @examples
#' toy_data <- simulate_csmart_data(seed = 11)
#' head(toy_data)
#' @export
simulate_csmart_data <- function(n_clusters = 24,
                                 cluster_size = 12,
                                 rerandomization = c(
                                   "nonresponder",
                                   "responder",
                                   "full"
                                 ),
                                 p_rerand = 0.6,
                                 icc = 0.1,
                                 seed = NULL) {
  rerandomization <- match.arg(rerandomization)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(cluster_size) == 1L) {
    cluster_size <- rep(as.integer(cluster_size), n_clusters)
  }
  if (length(cluster_size) != n_clusters) {
    stop("`cluster_size` must have length 1 or `n_clusters`.", call. = FALSE)
  }

  cluster_id <- seq_len(n_clusters)
  X1 <- sample(c(-1, 1), n_clusters, replace = TRUE)
  A1 <- sample(c(-1, 1), n_clusters, replace = TRUE)
  p_x2 <- plogis(0.3 * X1 + 0.2 * A1)
  X2 <- ifelse(stats::runif(n_clusters) < p_x2, 1, -1)

  if (rerandomization == "full") {
    rerandomized <- rep(1L, n_clusters)
    response_status <- rep(NA_integer_, n_clusters)
  } else {
    response_status <- stats::rbinom(n_clusters, size = 1, prob = 1 - p_rerand)
    rerandomized <- if (rerandomization == "nonresponder") {
      1L - response_status
    } else {
      response_status
    }
  }

  A2 <- rep(NA_real_, n_clusters)
  A2[rerandomized == 1L] <- sample(c(-1, 1), sum(rerandomized == 1L), replace = TRUE)
  A2_effective <- ifelse(is.na(A2), A1, A2)

  cluster_effect <- stats::rnorm(n_clusters, sd = sqrt(icc))

  dat_list <- vector("list", n_clusters)
  for (i in seq_len(n_clusters)) {
    n_i <- cluster_size[i]
    individual_error <- stats::rnorm(n_i, sd = sqrt(max(1 - icc, 1e-8)))
    mu <- 0.5 +
      0.6 * X1[i] +
      0.3 * A1[i] +
      0.15 * X1[i] * A1[i] +
      0.35 * A2_effective[i] +
      0.25 * A1[i] * A2_effective[i] +
      0.30 * X2[i] * A2_effective[i]

    dat_list[[i]] <- data.frame(
      cluster_id = rep(cluster_id[i], n_i),
      patient_id = seq_len(n_i),
      X1 = rep(X1[i], n_i),
      A1 = rep(A1[i], n_i),
      X2 = rep(X2[i], n_i),
      A2 = rep(A2[i], n_i),
      response_status = rep(response_status[i], n_i),
      rerandomized = rep(rerandomized[i], n_i),
      Y = rep(mu + cluster_effect[i], n_i) + individual_error
    )
  }

  do.call(rbind, dat_list)
}
