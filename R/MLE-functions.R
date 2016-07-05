kappahat = function(m, y, theta) {
  pfun = function(y, x)
    (y - x)*(y > x)

  mfun = function(y, x)
    (x - y)*(y < x)

  return(((1 / m * (sum(mfun(y, theta)))) / (1 / m * (sum(pfun(y, theta)))))^(1/4))
}

sigmahat1 = function(m, y, theta) {
  pfun = function(y, x)
    (y - x)*(y > x)

  mfun = function(y, x)
    (x - y)*(y < x)

  return(sqrt(2) * (1 / m * sum(pfun(y, theta)))^(1/4) *
           (1 / m * sum(mfun(y, theta))))^(1/4) *
    (sqrt(1 / m * sum(pfun(y, theta))) + sqrt(1 / m * sum(mfun(y, theta))))
}

sigmahat2 = function(m, y, theta, kappa) {
  pfun = function(y, x)
    (y - x)*(y > x)

  mfun = function(y, x)
    (x - y)*(y < x)

  return(sqrt(2) / m * (sum(pfun(y, theta)) * kappa + sum(mfun(y, theta)) / kappa))
}
