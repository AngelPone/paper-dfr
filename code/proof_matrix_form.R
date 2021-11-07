# This R code is used to check the validity of the matrix form


lambda = runif(1)

time_window = 100

# pi_hat: q * T
pi_hat = matrix(rnorm(q*time_window), nrow = q)
# A: r * q
A = matrix(rnorm(r*q), ncol = q)
distance = matrix(rnorm(r*q), ncol = q)
# Z: q * T
Z = matrix(0, nrow = r, ncol = time_window)
for (i in 1:time_window){
  Z[sample(1:r, 1), i] = 1
}

# scalar form
total = 0
for (i in 1:time_window){
  for (j in 1:r){
    total = total + (sum(A[j,] * pi_hat[, i]) - Z[j, i])^2
  }
}
total = total / time_window
total = total + lambda * sum(A * distance)


# Matrix form

vec_A <- matrix(t(A), ncol = 1)
Q <-  pi_hat %*% t(pi_hat)
Q <- bdiag(Q, Q, Q, Q)
vec_Z <- matrix(Z, nrow = 1)
C <- NULL

for (i in 1:time_window){
  C <- rbind(C, bdiag(replicate(r, matrix(pi_hat[,i], nrow = 1), simplify = FALSE)))
}
total_matrix = (t(vec_A) %*% Q/time_window) %*% vec_A - 
  2*(vec_Z %*% C/time_window) %*% vec_A + 
  lambda * matrix(t(distance), nrow = 1) %*% vec_A +
  vec_Z %*% t(vec_Z)/time_window


print((total - total_matrix[1][1]) < 1e-10)

total

total_matrix

debug(opt_fun)
opt_fun(pi_hat, Z, distance, lambda)
