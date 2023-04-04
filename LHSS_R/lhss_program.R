
set.seed(1)

N <- 1000
x_nu1 <- rnorm(N, 0, 1)
x_de1 <- rnorm(N, 0, 1)

x_nu2 <- rnorm(N, 0, 0.5)
x_de2 <- sample(rep(c(-1, 1), times = N/2)) + rnorm(N, 0, 1)

X_nu <- cbind(x_nu1, x_nu2)
X_de <- cbind(x_de1, x_de2)

p <- ncol(X_nu)
n_nu <- n_de <- N
m <- 1
nBasis <- 100

# X_nu <- rnorm(150) |> matrix(25)
# X_de <- rnorm(150) |> matrix(25)
# 
# p <- ncol(X_nu)
# n_nu <- nrow(X_nu)
# n_de <- nrow(X_de)
# 
# m <- 3

U  <- matrix(1, p, m)
QR <- pracma::householder(U)
Q  <- QR$Q * sign(QR$R[1])
U  <- Q[,1:m, drop = F]

centers <- X_nu[sample(n_nu, nBasis), ]

temp <- X_nu %*% U
sigma <- (densityratio:::distX(X_nu, nrow(X_nu), ncol(X_nu))/2) |> 
  median() |> 
  sqrt()

sigma <- 0.5

phi_nu <- densityratio:::kernel_gaussian(X_nu, centers, sigma)
phi_de <- densityratio:::kernel_gaussian(X_de, centers, sigma)
Hhat   <- crossprod(phi_de) / n_de
one    <- diag(0.1, ncol(Hhat))
hhat   <- colMeans(phi_nu)

alphat <- solve(Hhat + one) %*% hhat
alphah <- pmax(0, alphat)
PD_opt <- crossprod(hhat, alphah) - 0.5
U_opt  <- U

decPDcount <- 0
conv <- FALSE
maxit <- 200
iter <- 0

while(!conv) {
  
  iter <- iter+1
  cat(paste0("\r Iteration: ", iter))
  print(U)
  
  QR <- pracma::householder(U)
  Q  <- QR$Q * sign(QR$R[1])
  U  <- Q[, 1:m, drop = F]
  V  <- Q[, (m+1):p, drop = F]
  
  u_x_ce <- centers %*% U
  h <- rep(0, p)
  u_x_nu <- X_nu %*% U
  u_x_de <- X_de %*% U
  
  dPd1 <- matrix(0, p, m)
  Ktemp <- densityratio:::kernel_gaussian(u_x_nu, u_x_ce, sigma)
  
  for (i in 1:nBasis) {
    temp <- -(u_x_nu - rep(1, n_nu) %*% t(u_x_ce[i, ])) * (Ktemp[, i] %*% t(rep(1, m)))
    temp1 <- (X_nu - rep(1, n_nu) %*% t(centers[i, ]))
    h <- crossprod(temp, temp1)
    
    dPd1 <- dPd1 + t(h) * alphah[i]
  }
  
  dPd1 <- dPd1 / n_nu / sigma
  
  Ktemp <- densityratio:::kernel_gaussian(u_x_de, u_x_ce, sigma)
  dPd2 <- matrix(0, p, m)
  
  for (i in 1:nBasis) {
    temp11 <- (u_x_de - rep(1, n_de) %*% t(u_x_ce[i, ])) * (Ktemp[, i] %*% t(rep(1, m)))
    tempx1 <- (X_de - rep(1, n_de) %*% t(centers[i, ]))
    
    for (j in 1:m) {
      T11 <- temp11[,j] %*% t(rep(1, nBasis)) * Ktemp
      dPd2[, j] <- dPd2[, j] - t(tempx1) %*% T11 %*% alphah * alphah[i] * 2
    }
  }
  
  dPd2 <- dPd2 / n_de / sigma
  
  dPd <- t(dPd1 - dPd2/2)
  
  dM <- rbind(
    cbind(matrix(0, m, m), dPd %*% V),
    cbind(-t(dPd %*% V), matrix(0, p - m, p - m))
  )
  
  M <- dM * 1/maxit * (maxit - (iter - 1))
  eM <- expm::expm(M)
  U <- cbind(diag(m), matrix(0, m, p - m)) %*% eM %*% rbind(t(U), t(V))
  U <- t(U)
  
  temp <- X_nu %*% U
  med <- (densityratio:::distX(temp, nrow(temp), ncol(temp))/2) |> 
    median() |> 
    sqrt()
  
  phi_nu <- densityratio:::kernel_gaussian(X_nu %*% U, centers %*% U, med)
  phi_de <- densityratio:::kernel_gaussian(X_de %*% U, centers %*% U, med)
  Hhat   <- crossprod(phi_de) / n_de
  one    <- diag(1, ncol(Hhat))
  hhat   <- colMeans(phi_nu)
  
  alphat <- solve(Hhat + one) %*% hhat
  alphah <- pmax(0, alphat)
  PD <- crossprod(hhat, alphah) - 0.5
  
  if (PD > PD_opt) {
    U_opt <- U
    alphah_opt <- alphah
    PD_opt <- PD
    sigmatopt <- med
    decPDcount <- max(0, decPDcount - 1)
  } else {
    decPDcount <- decPDcount + 1
  }
  
  if (iter == maxit | decPDcount > 20) {
    conv <- T
  }
}

alphah

library(ggplot2)
library(dplyr)

bind_rows(nu = X_nu |> data.frame() |> select(V = x_nu1, U = x_nu2),
          de = X_de |> data.frame() |> select(V = x_de1, U = x_de2),
          .id = "Samples") |>
  ggplot(aes(x = U, y = V, col = Samples)) +
  geom_point() +
  theme_minimal() +
  scale_color_viridis_d(alpha = .5)


xmat <- tidyr::expand_grid(x1 = seq(-5, 5, by = 0.05),
                           x2 = seq(-5, 5, by = 0.05)) |>
  as.matrix()

K_test <- densityratio:::kernel_gaussian(xmat %*% U, y = centers %*% U, sigma = sigma)

xmat |>
  data.frame() |>
  mutate(r = densityratio:::kernel_gaussian(xmat, centers, sigma) %*% 
           densityratio::ulsif(X_nu, X_de, sigma, 1, centers = centers)) |>
  ggplot(aes(x = x2, y = x1, z = r)) +
  stat_contour_filled()

xmat |>
  data.frame() |>
  mutate(r = K_test %*% alphah) |>
  ggplot(aes(x = x2, y = x1, z = r)) +
  stat_contour_filled()
