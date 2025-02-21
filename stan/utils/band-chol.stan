matrix band_chol(matrix Q, int bw) {
  int N = rows(Q);
  matrix[N, N] L = rep_matrix(0, N, N);
  for (j in 1:N) {
    int lambda = min(j + bw, N);
    vector[N] v = rep_vector(0, N);
    v[j:lambda] = Q[j:lambda, j];
    int k_start = max(1, j - bw);
    int k_end = j - 1;
    for (k in k_start:k_end) {
      int i = min(k + bw, N);
      v[j:i] = v[j:i] - L[j:i, k] * L[j, k];
    }
    L[j:lambda, j] = v[j:lambda] * inv_sqrt(v[j]);
  }
  return L;
}
vector backsolve_bchol(matrix L, vector z, int bw) {
  int N = rows(L);
  vector[N] x = z;
  for (i in 1:N) {
    int k_start = max(1, i - bw);
    int k_end = i - 1;
    for (j in k_start:k_end) {
      x[i] = x[i] - L[i, j] * x[j];
    }
    x[i] = x[i] / L[i, i];
  }
  return x;
}
matrix solve_bchol(matrix L, int bw) {
  int N = rows(L);
  matrix[N, N] Linv;
  for (j in 1:N) {
    vector[N] aux = rep_vector(0.0, N);
    aux[j] = 1.0;
    Linv[1:, j] = backsolve_bchol(L, aux, bw);
  }
  return Linv;
}

matrix ar1_prec(real rho, int size) {
  matrix[size, size] Q = rep_matrix(0.0, size, size);
  for (i in 1:size) {
    int k_start = max(1, i - 1);
    int k_end = min(size, i + 1);
    for (j in k_start:k_end) {
      if (j == 1 || j == size) {
        Q[i, j] = i == j ? 1 : - rho;
      } else {
        Q[i, j] = i == j ? 1 + square(rho) : - rho;
      }
    }
  }
  return Q;
}
