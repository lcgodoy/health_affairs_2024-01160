real pcp_ar0_lpdf(real x, real alpha, real u) {
  real aux = - log(1 - square(x));
  real theta = - log(alpha) / sqrt(- log(1 - square(u)));
  return log(theta) - log(2) - theta * sqrt(aux) +
    log(abs(x)) + aux - 0.5 * log(aux);
}

// for the future
// real pcp_ar1_lpdf(real x, real alpha, real u) {
// }
