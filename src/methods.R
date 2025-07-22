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
fit_gamma <- function(obs, obs_thresh = 0.85, method = "mle") {
  obs <- obs[obs > obs_thresh] - obs_thresh
  obs <- as.numeric(na.omit(obs))
  if(length(obs) > 5) fitdistrplus::fitdist(obs, distr = "gamma", method = method)
  else NULL
}

# Weather generator

