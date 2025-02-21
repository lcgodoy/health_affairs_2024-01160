##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Auxilliary function for prior sensitivity
##' @param x a list subset of a `powerscaled_sequence` object
##' @param fun a `function` having a `vector` as an input and returning a
##'   `vector` with summary statistics
##' @param cnames a `character` vector of the same size as the output of `fun`
##' @param tp a `character`
##' @return a `tibble`.
##' @author Lucas Godoy
aux_sens <- function(x, fun, cnames = NULL, tp) {
  if ("powerscaling" %in% names(x)) {
    alpha <- x$powerscaling$alpha
    out <- apply(x[["draws"]], 2, fun)
  } else {
    alpha <- 1
    out <- apply(x, 2, fun)
  }
  if (NCOL(out) > 1) {
    out <- t(out)
    par <- rownames(out)
  } else {
    par <- names(out)
    out <- matrix(out, ncol = 1)
  }
  if (!is.null(cnames))
    colnames(out) <- cnames
  out <- dplyr::as_tibble(out) |>
    dplyr::mutate(alpha = alpha,
                  type = tp,
                  parameter = par,
                  .before = 1)
  out <- out |>
    dplyr::filter(!grepl("^\\.", out[["parameter"]]))
  rownames(out) <- NULL
  return(out)
}

##' @title Functions for prior-sensitivity analysis
##' @param pwr_sq a `powerscaled_sequence` object
##' @param stats a function that summarises `mcmc` samples
##' @param nms names of the outputs of `stats`
##' @return a `tibble`
##' @author Lucas Godoy
sens_fun <- function(pwr_sq, stats, nms = NULL) {
  prior_ests <- lapply(pwr_sq$prior_scaled[[1]],
                       aux_sens,
                       fun = stats,
                       cnames = nms,
                       tp = "prior")
  lik_ests <- lapply(pwr_sq$likelihood_scaled[[1]],
                     aux_sens,
                     fun = stats,
                     cnames = nms,
                     tp = "likelihood")
  observed <- aux_sens(pwr_sq$base_draws,
                       fun = stats,
                       cnames = nms,
                       tp = "observed")
  lik_ests <-
    lik_ests |>
    dplyr::bind_rows(
               dplyr::mutate(observed, type = "likelihood")
           )
  prior_ests <-
    prior_ests |>
    dplyr::bind_rows(
               dplyr::mutate(observed, type = "prior")
           )
  return(dplyr::bind_rows(lik_ests, prior_ests))
}
