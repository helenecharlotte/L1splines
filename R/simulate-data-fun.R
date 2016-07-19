##--- generate sample times
#' @title Generate sample times.
#'
#' Function to generate different kinds of sample times.
#'
#' @param type Numeric value (1,2).
#' @param m Length of output vector.
#' @param seed If certain seed wanted.
#' @examples
#' generate.t()
#' @export
generate.t = function(type = 1, m = 15, seed = sample(2500, 1)) {
  set.seed(seed)
  if (type == 1)
    t = seq(0, 1, length = m)
  if (type == 2)
    t = c(0, sort(runif(m-2)), 1)
  return(t)
}

##--- generate observations of function to be fitted
#' @title Generate observations (underlying functions).
#'
#' Function to generate different kinds of underlying functions.
#'
#' @param type Numeric value (1,2,...,10).
#' @param seed If certain seed wanted.
#' @examples
#' y = generate.y(); plot(y)
#' y = generate.y(type = 2); plot(y)
#' y = generate.y(type = 3); plot(y)
#' y = generate.y(type = 4); plot(y)
#' @export
generate.y = function(type = 1, seed = sample(2500, 1)) {
  set.seed(seed)
  a = runif(1); b = runif(1)
  if (type == 1)
    y = function(t)
      cos(8 * t^2) - t^3
  if (type == 2) ## monotonely increasing
    y = function(t)
      splinefun(c(0, a, 1),
                c(0, b, 1), method = "hyman")(t)
  if (type == 3) ## convex
    y = function(t)
      (t - a)^2
  if (type == 4) {
    n = sample(10,1) + 45
    b = sort(runif(2))
    a = runif(n, min = b[1], max = b[2])
    tpoints = seq(0, 1, length = 100)
    thet    = sapply(1:n, function(i)
      dnorm(tpoints, mean = a[i], sd = 0.01))
    ypoints = sapply(1:length(tpoints), function(i)
      sum(thet[i, ]))
    y = function(t) approxfun(tpoints, ypoints, rule = 2)(t)
  }
  if (type == 5) {
    n = sample(10,1) + 45
    b = sort(runif(2))
    a = runif(n, min = b[1], max = b[2])
    tpoints = seq(0, 1, length = 100)
    thet    = sapply(1:n, function(i)
      dnorm(tpoints, mean = a[i], sd = 0.01))
    ypoints = c(rep(0.3, 45), rep(1, 10), rep(0.3, 20),
                rep(1, 15), rep(0.3, 10)) * sapply(1:length(tpoints), function(i)
      sum(thet[i, ]))
    y = function(t) approxfun(tpoints, ypoints, rule = 2)(t)
  }
  if (type == 6)
    y = function(t)
      0.01 * dnorm(t, mean = 0.2, sd = 0.3) +
    0.0025 * dnorm(t, mean = max(a, 0.6), sd = 0.05) +
    0.0025 * dnorm(t, mean = max(b, 0.5), sd = 0.05)
  if (type == 7)
    y = function(t)
      0.8 * dnorm(t, mean = 0.2, sd = 0.3)
  if (type == 8)
    y = function(t)
      log(t + 1.0)
  if (type == 9)
    y = function(t)
      exp(-abs(t-a)/(b/38))/(2*b)
  x = seq(0, 1, length = 100)
  if (type == 10)
    y = approxfun(x,  (abs(x - 1/2) > 0.001) * 1 / (15 * (x - 1/2)))
  if (type == 11) ## convex again
    y = function(t)
      t^2 - t
  if (type == 12)
    y = function(t)
      sin(t*30)
  if (type == 13)
    y = function(t)
      1 - abs(t-2/3) - abs(t-1/3)
  if (type == 14)
    y = function(t) {
      i = findInterval(t, c(0, a, 1), all.inside = TRUE)
      return(c(a, b)[i])
    }
  ynorm = sapply(seq(0, 1, length = 100), function(t) y(t) - y(0))
  return(function(t) 2.5 * (y(t) - y(0)) / (max(ynorm) - min(ynorm)))
}


##--- generate noise
#' @title Generate random noise.
#'
#' Function to generate different kinds of random noise.
#'
#' @param t Numeric vector where to evaluate noise.
#' @param type Numeric value (1,...,8).
#' @param mu Numeric value.
#' @param theta Numeric value.
#' @param kappa Numeric value.
#' @param sigma Numeric value.
#' @param sigma2 Numeric value. Smaller SE for contaminated noise.
#' @param seed If certain seed wanted.
#' @param tfun For generating noise with dynamically changing parameters.
#' @examples
#' t = generate.t(m = 50)
#' w = generate.noise(t); plot(t, w , type="h", main = "Gaussian noise")
#' w = generate.noise(t, type = 2); plot(t, w , type="h", main = "Laplace noise")
#' w = generate.noise(t, type = 3); plot(t, w , type="h", main = "Contaminated noise")
#' w = generate.noise(t, type = 4); plot(t, w , type="h", main = "Uniform noise")
#' @export
generate.noise = function(t, type = 1, mu = 0, theta = 0, kappa = Inf,
                          sigma = 0.2, sigma2 = 1e-4,
                          seed = sample(2500, 1),
                          tfun = function(t) t / 4) {
  set.seed(seed)
  m = length(t)
  if (type == 1)
    noise = rnorm(m, mean = mu, sd = sigma)
  if (type == 2)
    noise = rLaplace(m, mu = mu, theta = theta, kappa = kappa, sigma = sigma)
  if (type == 3) {
    w     = rbinom(m, 1, 0.35);
    noise = w * rnorm(m, sd = sigma) + (1-w) * rnorm(m, sd = sigma2)
  }
  if (type == 4)
    noise = runif(m, min = -sigma, max = sigma)
  if (type == 5)
    noise = c(1, rep(0, length(t)-2), 1) * runif(m, min = -sigma, max = sigma) +
      rnorm(m, sd = sigma / 2)
  if (type == 6) {
    noise = sapply(t, function(t)
      rLaplace(1, mu = mu, theta = theta, kappa = kappa, sigma = sigma * tfun(t + 0.01)))
    attr(noise, "sigma") = sigma * tfun(t + 0.01)
  }
  if (type == 7) {
    noise = sapply(t, function(t)
      rLaplace(1, theta = 0, kappa = kappa * tfun(t + 0.01), sigma = sigma))
    attr(noise, "kappa") = kappa * tfun(t + 0.01)
  }

  return(noise)
}

