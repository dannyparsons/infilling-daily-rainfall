# Direct replacement

direct_replace <- function(obs, est) {
  ix <- is.na(obs)
  obs[ix] <- est[ix]
  obs
}

# Threshold adjustment
thresh_calc <- function(obs, est, obs_thresh = 0.85) {
  id <- !is.na(obs) & !is.na(est)
  p_obs <- mean(obs[id] > obs_thresh, na.rm = TRUE)
  quantile(est[id], 1 - p_obs)[[1]]
}

# Local Intensity Scaling

# Quantile Mapping

# Weather generator
