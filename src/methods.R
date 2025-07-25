# Direct replacement

direct_replace <- function(obs, est) {
  ix <- is.na(obs)
  obs[ix] <- est[ix]
  obs
}

# Threshold adjustment
thresh <- function(obs, est, obs_thresh = 0.85) {
  id <- !is.na(obs) & !is.na(est)
  p_obs <- mean(obs[id] > obs_thresh, na.rm = TRUE)
  quantile(est[id], 1 - p_obs)[[1]]
}

s_loci <- function(obs, est, obs_thresh = 0.85, est_thresh = NULL) {
  if (is.null(est_thresh)) est_thresh <- thresh(obs, est, obs_thresh)
  id <- !is.na(obs) & !is.na(est)
  s_obs <- mean(obs[id & obs > obs_thresh]) - obs_thresh
  s_est <- mean(est[id & est > est_thresh]) - est_thresh
  s <- s_obs / s_est
  s
}

# Local Intensity Scaling
loci <- function(obs, est, obs_thresh = 0.85) {
  est_thresh <- thresh(obs, est, obs_thresh)
  s <- s_loci(obs, est, obs_thresh, est_thresh)
  est_loci <- obs_thresh + s * (est - est_thresh)
  est_loci <- pmax(est_loci, 0)
  est_loci
}

# Quantile Mapping
fit_gamma <- function(obs, obs_thresh = 0.85, 
                      method = c("mle", "mme", "qme", "mge", "mse")) {
  obs <- obs[obs > obs_thresh] - obs_thresh
  obs <- as.numeric(na.omit(obs))
  if(length(obs) > 5) fitdistrplus::fitdist(obs, distr = "gamma", method = method)
  else NULL
}

fit_empirical <- function(obs, obs_thresh = 0.85) {
  obs <- obs[obs > obs_thresh] - obs_thresh
  ecdf(obs)
}

quantile_mapping <- function(obs, est, obs_thresh = 0.85, 
                             qm_method = c("gamma", "empirical"), 
                             gamma_method = c("mle", "mme", "qme", 
                                              "mge", "mse")) {
  est_thresh <- thresh(obs, est, obs_thresh)
  qm_method <- match.arg(qm_method)
  if (qm_method == "gamma") {
    obs_gamma <- fit_gamma(obs, obs_thresh, method = gamma_method)
    est_gamma <- fit_gamma(est, est_thresh, method = gamma_method)
    obs_shape = obs_gamma[["estimate"]][["shape"]]
    obs_rate = obs_gamma[["estimate"]][["rate"]]
    est_shape = est_gamma[["estimate"]][["shape"]]
    est_rate = est_gamma[["estimate"]][["rate"]]
    est_qm <- if_else(est <= est_thresh, 0,
                      qgamma(pgamma(est - est_thresh, shape = est_shape, 
                                    rate = est_rate),
                             shape = obs_shape, rate = obs_rate) + obs_thresh)
  } else if (qm_method == "empirical") {
    obs_rain <- obs[obs > obs_thresh] - obs_thresh
    est_ecdf <- fit_empirical(est, est_thresh)
    est_qm <- if_else(est <= est_thresh, 0,
                      as.numeric(quantile(obs_rain, est_ecdf(est - est_thresh), 
                               na.rm = TRUE)) + obs_thresh)
  }
  return(est_qm)
}

# Weather generator

