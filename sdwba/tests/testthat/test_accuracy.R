context("Accuracy")

sound.speed <- 1456

generic.krill.m <- generic.krill.McGeehee1998
generic.krill.m[c("x", "y", "z", "a")] <- generic.krill.m[c("x", "y", "z", "a")] * 1e-3
krill <- with(generic.krill.m, Scatterer(x, y, z, a, g, h))


freqs <- c(38, 120, 200, 333) * 1e3

for (f in freqs) {
  k <- c(0, 0, -f) * 2 * pi / sound.speed
  print(paste(f/1000, target.strength(krill, k)))
}

freq.spec <- frequency.spectrum(krill, 12e3, 400e3, sound.speed, nfreq=120)
plot(freq.spec$freq, freq.spec$TS, ty='l')

angle.spec <- tilt.spectrum(krill, -90, 270, 120e3, sound.speed, nangle=180)
angle.spec.c <- tilt.spectrum(krill, -90, 270, 120e3, sound.speed, nangle=180,
                            method="continuous")
plot(angle.spec$angle, angle.spec$TS, ty='l')
lines(angle.spec.c$angle, angle.spec.c$TS, col='blue')
max(angle.spec$TS)



angles <- -90:270
sigma <- rep(0, length(angles))
k <- c(0, 0, -120e3) * 2 * pi / sound.speed
phase.sd <- sqrt(2) / 2

for (i in 1:length(angles)) {
  print(i)
  sigma[i] <- backscatter.xsection.ensemble(rotate(krill, angles[i]), k, phase.sd)$mean
}
TS <- 10 * log10(sigma)
plot(angles, TS, ty='l')




## Fluid-sphere

besselj <- function(x, nu) sqrt(0.5 * pi / x) * besselJ(x, nu + 0.5)
bessely <- function(x, nu) sqrt(0.5 * pi / x) * besselY(x, nu + 0.5)
alpha <- function(x, m) m * besselj(x, m-1) - (m + 1) * besselj(x, m+1)
beta <- function(x, m) m * bessely(x, m-1) - (m - 1) * bessely(x, m+1)


backscatter.sphere <- function(a, k, g, h, M=10) {
  s <- 0
  k1 <- k / h
  for (m in 0:M) {
    C <- alpha(k1*a, m) / alpha(k*a, m) * bessely(k*a, m) / besselj(k1*a, m) -
      beta(k*a, m) / alpha(k*a, m) * g * h
    C <- C / (alpha(k1*a, m) / alpha(k*a, m) * besselj(k*a, m) / besselj(k1*a, m) - g * h)
    s <- s + (-1)^m * (2 * m + 1) / (1 + 1i * C)
  }
  R <- 2 / (k * a) * abs(s)
  return(R^2 * (2 * pi * a^2))
}

g <- 1.0
h <- 1.1
a <- 1
kk <- seq(0.4, 4, length=30) / a
hh <- seq(0.5, 2, length=30)
R <- matrix(0, length(kk), length(hh))
for (i in 1:length(kk)) {
  for (j in 1:length(hh)) {
    R[i, j] <- backscatter.sphere(a, kk[i], g, hh[j]) / (2 * pi * a^2)
  }
}
persp(kk, hh, R, theta=130)





