library(numDeriv)

# convenience functions for vectors and trig

#' Calculate the norm (magnitude) of a vector.
#'
#' @param x A vector
#' @return sqrt(sum(x^2))
#' @examples
#' x <- c(3, 4)
#' norm(x)
#' [1] 5
norm <- function(x) sqrt(sum(x^2))

#' Calculate the dot product of two vectors.
#'
#' @param x,y Two vectors with the same length.
#' @return The dot (inner) product of x and y
dot <- function(x, y) sum(x * y)

#' Convert between degrees and radians.
#'
#' @param x Angle to convert
#' @usage
#' rad2deg <- function(x) x * 180 / pi
#' deg2rad <- function(x) x * pi / 180
#' @examples
#' rad2deg(pi)
#' [1] 180
#' deg2rad(180)
#' [1] 3.141593
rad2deg <- function(x) x * 180 / pi
deg2rad <- function(x) x * pi / 180


#' Scatterer Object
#'
#' @description
#' This function creates Scatterer objects, which encapsulate the shape and
#' material properties of a fluid-like object with a distorted cylindrical shape
#' which scatters sound.  The shape is approximated by a series of short cylinders,
#' which can have varying radii, densities, and sound speeds.
#' @param x,y,z Vectors defining the scatterer's centerline.  x and y are horizontal, z is up.
#' @param a Vector of radii
#' @param g Vector of density contrasts
#' @param h Vector of sound-speed contrasts
Scatterer <- function(x, y, z, a, g, h) {
  scatterer <- data.frame(x=x, y=y, z=z, a=a, g=g, h=h)
  class(scatterer) <- append(class(scatterer), "Scatterer")
  return(scatterer)
}

#' Load a Scatterer from a CSV file.
#'
#' @description
#' Reads in a scatterer from a CSV file with columns matching the arguments to
#' Scatterer().
#' @param filename Path to a CSV file. This file must have columns named "x", "y",
#' "z", "a", "g", and "h" defining the scatterer's shape and material properties.
#' @param ... Additional arguments passed to read.csv
#' @details This function is a wrapper around read.csv, allowing a Scatterer to
#' be loaded directly from a CSV file.  The file must contain columns with the
#' headers "x", "y", "z", "a", "g", and "h".  It may contain other columns as well,
#' but these will be discarded.
#' @return A Scatterer object.
load.scatterer <- function(filename, ...) {
  df <- read.csv(filename, ...)
  return(Scatterer(x=df$x, y=df$y, z=df$z, a=df$a, g=df$g, h=df$h))
}


save.scatterer <- function(scatterer, ...)  UseMethod("save.scatterer", scatterer)

#' Save a Scatterer object to file.
#'
#' @description Writes a Scatterer object to a CSV file.
#' @param scatterer The Scatterer object to save
#' @param filename Filename to save to
#' @param ... Additional arguments passed to write.csv
save.scatterer.Scatterer <- function(scatterer, filename, ...) {
  write.csv(scatterer, filename, ...)
}

rotate <- function(scatterer, ...) UseMethod("rotate", scatterer)

#' Rotate a Scatterer in space.
#'
#' @description Tilts and/or rotates a scatterer, returning a new scatterer object.
#' @param scatterer A Scatterer object
#' @param roll Angle to roll the scatterer, in degrees. Defaults to 0.
#' @param tilt Angle to tilt the scatterer, in degrees. Defaults to 0.
#' @param yaw Angle to yaw the scatterer, in degrees. Defaults to 0.
#' @return A Scatterer with the same shape and properties, but a new orientation.
#' @details
#' The roll, tilt, and yaw refer to rotations around the x, y, and z axes,
#' respectively. They are applied in that order.
rotate.Scatterer <- function(scatterer, roll=0, tilt=0, yaw=0) {
  tilt <- deg2rad(tilt)
  roll <- deg2rad(roll)
  yaw <- deg2rad(yaw)
  Rx <- matrix(c(1, 0, 0, 0, cos(roll), sin(roll), 0, -sin(roll), cos(roll)), nrow=3)
  Ry <- matrix(c(cos(tilt), 0, -sin(tilt), 0, 1, 0, sin(tilt), 0, cos(tilt)), nrow=3)
  Rz <- matrix(c(cos(yaw), sin(yaw), 0,  -sin(yaw), cos(yaw), 0, 0, 0, 1), nrow=3)
  R <- Rz %*% Ry %*% Rx
  for (i in 1:nrow(scatterer)) {
    scatterer[i, c('x', 'y', 'z')] <- R %*% unlist(scatterer[i, c('x', 'y', 'z')])
  }
  return(scatterer)
}

rescale <- function(scatterer, ...) UseMethod("rescale", scatterer)
#' Resize a Scatterer
#'
#' This function stretches or resizes a scatterer.
#' @param scatterer A Scatterer object
#' @param scale Overall scaling factor to apply.  Each dimension of the scatterer
#' will be multiplied by this number.  Defaults to 1.
#' @param radius Factor by which to multiply the Scatterer's radii.  Defaults to 1.
#' @param x,y,z Individual scaling factors to apply to the x, y, and z coordinates
#' of the Scatterer's centerline.  All default to 1.
#' @details All arguments to this function are applied; all default to 1 (i.e.,
#' no stretching). You should only specify the ones you wish to change.
#' @return A rescaled Scatterer object.
rescale.Scatterer <- function(scatterer, scale=1, radius=1, x=1, y=1, z=1) {
  scatterer$x <- scatterer$x * scale * x
  scatterer$y <- scatterer$y * scale * y
  scatterer$z <- scatterer$z * scale * z
  scatterer$a <- scatterer$a * scale * radius
  return(scatterer)
}

total.length <- function(scatterer, ...) UseMethod("total.length", scatterer)
total.length.Scatterer <- function(scatterer) {
  n <- nrow(scatterer)
  return(norm(scatterer[1, c("x", "y", "z")] - scatterer[n, c("x", "y", "z")]))
}

resample <- function(scatterer, ...) UseMethod("resample", scatterer)

resample.Scatterer <- function(scatterer, f, f0) {
  N0 <- nrow(scatterer)
  L0 <- total.length(scatterer)
  N <- N0 * f / f0

  interps <- interpolation.functions(scatterer)
  attach(interps)

  detach(interps)

}


plot <- function(scatterer, ...) UseMethod("plot", scatterer)

#' Visualize the shape of a Scatterer.
#'
#' @param scatterer A scatterer object
#'
plot.Scatterer <- function(scatterer) {
  xlim <- range(scatterer$x)
  ylim <- range(scatterer$z)
  w <- diff(xlim)
  h <- diff(ylim)
  xlim <- xlim + 0.1 * w * c(-1, 1)
  ylim <- ylim + 0.1 * h * c(-1, 1)
  plot.new()
  plot.window(xlim, ylim, asp=1, xaxs='r', yaxs='r')
  lines(scatterer$x, scatterer$z, ty="b")

  for (i in 1:(nrow(scatterer) - 1)) {
    dx <- scatterer$x[i+1] - scatterer$x[i]
    dy <- scatterer$z[i+1] - scatterer$z[i]
    theta <- atan2(dx, dy)
    x.top.1 <- scatterer$x[i] - scatterer$a[i] * cos(theta)
    x.bot.1 <- scatterer$x[i] + scatterer$a[i] * cos(theta)
    y.top.1 <- scatterer$z[i] + scatterer$a[i] * sin(theta)
    y.bot.1 <- scatterer$z[i] - scatterer$a[i] * sin(theta)

    x.top.2 <- scatterer$x[i + 1] - scatterer$a[i + 1] * cos(theta)
    x.bot.2 <- scatterer$x[i + 1] + scatterer$a[i + 1] * cos(theta)
    y.top.2 <- scatterer$z[i + 1] + scatterer$a[i + 1] * sin(theta)
    y.bot.2 <- scatterer$z[i + 1] - scatterer$a[i + 1] * sin(theta)

    lines(c(x.bot.1, x.top.1), c(y.bot.1, y.top.1), lty=3)
    lines(c(x.bot.2, x.top.2), c(y.bot.2, y.top.2), lty=3)
    lines(c(x.top.1, x.top.2), c(y.top.1, y.top.2))
    lines(c(x.bot.1, x.bot.2), c(y.bot.1, y.bot.2))
    axis(1)
    axis(2)
  }
}


form.function <- function(scatterer, k, phase.sd=0) {
  n <- nrow(scatterer)
  f.bs <- 0 + 0i
  for (j in 1:(n - 1)) {
    r1 <- c(scatterer$x[j], scatterer$y[j], scatterer$z[j])
    r2 <- c(scatterer$x[j + 1], scatterer$y[j + 1], scatterer$z[j + 1])
    alphatilt <- acos(dot(k, (r2 - r1)) / (norm(k) * norm(r2 - r1)))
    betatilt <- abs(alphatilt - pi/2)

    # Define the function to integrate over this segment
    integrand <- function(s, real=TRUE) {
      rx <- s * (r2[1] - r1[1]) + r1[1]
      ry <- s * (r2[2] - r1[2]) + r1[2]
      rz <- s * (r2[3] - r1[3]) + r1[3]
      r <- c(rx, ry, rz)
      a <- s * (scatterer$a[j + 1] - scatterer$a[j]) + scatterer$a[j]
      h <- s * (scatterer$h[j + 1] - scatterer$h[j]) + scatterer$h[j]
      g <- s * (scatterer$g[j + 1] - scatterer$g[j]) + scatterer$g[j]
      gamgam <- 1 / (g * h^2) + 1 / g - 2
      k2 <- norm(k) / h

      if (abs(betatilt) < 1e-10) {
        bessy <- k2 * a / h
      } else {
        arg <- 2 * k2 * a / h * cos(betatilt)
        bessy <- besselJ(arg, 1) / cos(betatilt)
      }
      result <- norm(k) / 4 * gamgam * a * exp(2i * dot(k, r) / h) *
        bessy * norm(r2 - r1)
      if (real) {
        return(Re(result))
      } else {
        return(Im(result))
      }
    }
    integrand <- Vectorize(integrand)

    real.part <- integrate(integrand, 0, 1, real=TRUE)$value
    imag.part <- integrate(integrand, 0, 1, real=FALSE)$value
    segment.f.bs <- (real.part + imag.part) * exp(1i * rnorm(1, 0, deg2rad(phase.sd)))
    f.bs <- f.bs + segment.f.bs
  }
  return(f.bs)
}

interpolation.functions <- function(scatterer) {
  require(numDeriv)
  dx <- diff(scatterer$x)
  dy <- diff(scatterer$y)
  dz <- diff(scatterer$z)
  ss <- c(0, cumsum(sqrt(dx^2 + dy^2 + dz^2)))

  # Defining continuous interpolation functions
  x_fun <- splinefun(ss, scatterer$x)
  y_fun <- splinefun(ss, scatterer$y)
  z_fun <- splinefun(ss, scatterer$z)
  a_fun <- splinefun(ss, scatterer$a)
  g_fun <- splinefun(ss, scatterer$g)
  h_fun <- splinefun(ss, scatterer$h)
  r_fun <- function(s) c(x_fun(s), y_fun(s), z_fun(s))
  local_tangent <- function(s) jacobian(r_fun, s)

  functions <- list(x_fun=x_fun, y_fun=y_fun, z_fun=z_fun, a_fun=a_fun, g_fun=g_fun,
                    h_fun=h_fun, r_fun=r_fun, local_tangent=local_tangent)
  return(functions)
}

form.function.continuous <- function(scatterer, k, phase.sd=0) {
  require(numDeriv)
  # Note: phase.sd has a different meaning here than in the
  dx <- diff(scatterer$x)
  dy <- diff(scatterer$y)
  dz <- diff(scatterer$z)
  ss <- c(0, cumsum(sqrt(dx^2 + dy^2 + dz^2)))

  # Defining continuous interpolation functions
  x_fun <- splinefun(ss, scatterer$x)
  y_fun <- splinefun(ss, scatterer$y)
  z_fun <- splinefun(ss, scatterer$z)
  a_fun <- splinefun(ss, scatterer$a)
  g_fun <- splinefun(ss, scatterer$g)
  h_fun <- splinefun(ss, scatterer$h)
  r_fun <- function(s) c(x_fun(s), y_fun(s), z_fun(s))
  local_tangent <- function(s) jacobian(r_fun, s)

  # function to integrate along the length of the animal
  integrand <- function(s, k) {
    loc_tan <- local_tangent(s)
    beta <- acos(dot(k, loc_tan) / (norm(k) * norm(loc_tan)))
    beta <- abs(beta - pi/2)
    gamgam <- 1 / (g_fun(s) * h_fun(s)^2) + 1 / g_fun(s) - 2
    k2 <- norm(k) / g_fun(s)
    a <- a_fun(s)

    return(norm(k) / 4 * gamgam * exp(2i * dot(k, r_fun(s))) *
      a * besselJ(2 * k2 * a * cos(beta), 1) / cos(beta))
  }
  integrand.real <- Vectorize(function(s, k) Re(integrand(s, k)), vectorize.args=c("s"))
  integrand.imag <- Vectorize(function(s, k) Im(integrand(s, k)), vectorize.args=c("s"))

  f.bs.real <- integrate(integrand.real, 0, max(ss), k=k)$value
  f.bs.imag <- integrate(integrand.imag, 0, max(ss), k=k)$value
  f.bs <- (f.bs.real + f.bs.imag) #* exp(1i * rnorm(1, 0, phase.sd))
  return(f.bs)
}

#' Backscattering cross-section of a Scatterer.
#'
#' @description Calculates the backscattering cross-section of a scatterer at an arbitrary
#' frequency and incident angle.
#' @param scatterer A Scatterer object.
#' @param k A three-element vector, giving the x, y, and z components of the wavenumber
#' vector corresponding to the incident sound field.
#' @param phase.sd The standard deviation of the random phase component added to the form function
#' of each of the scatterer's segments.  In degrees, defaults to zero. Ignored if method is "continuous".
#' @param method The method to use when evaluating the line integral, defaults to "discrete".
#'
#' @details The two required arguments are scatterer and k.  The wavenumber vector has a magnitude
#' equal to the acoustic wavenumber, given by 2 * pi * f / c, where f is frequency and c is sound speed.
#' For a normal downward-looking echosounder, this is c(0, 0, -k).
#'
#' The two available integration methods are "discrete" and "continuous."  The discrete method
#' corresponds to that used in papers by McGeehee et al. (1998), and the series by Demer and Conti.
#' The continuous method is an experimental method which, instead of using the discretized outline of the
#' animal to define a series of short cylinders, uses it to interpolate the positions, radii, and
#' material properties along the animal's length, and integrate these functions continuously.
#' @return The backscattering cross section of the Scatterer.
backscatter.xsection <- function(scatterer, k, phase.sd=0, method=c("discrete", "continuous")) {
  method <- match.arg(method)
  if (method == "discrete") {
    f.bs <- form.function(scatterer, k, phase.sd)
  } else if (method == "continuous") {
    f.bs <- form.function.continuous(scatterer, k, phase.sd)
  } else {
    stop("method must be either 'discrete' or 'continuous'")
  }
  return(abs(f.bs)^2)
}

#' Target Strength of a Scatterer.
#'
#' This function is a simple wrapper around backscatter.xsection.
#' @param scatterer A Scatterer object.
#' @param k A three-element vector, giving the x, y, and z components of the wavenumber
#' vector corresponding to the incident sound field.
#' @param phase.sd The standard deviation of the random phase component added to the form function
#' of each of the scatterer's segments.  In degrees, defaults to zero. Ignored if method is "continuous".
#' @param method The method to use when evaluating the line integral, defaults to "discrete".
#'
#' @return The target strength (i.e., 10 * log10(backscatter.xsection(scatterer, k)))
#' @seealso backscatter.xsection
target.strength <- function(scatterer, k, phase.sd=0,
               method=c("discrete", "continuous")) {
  sigma.bs <- backscatter.xsection(scatterer, k, phase.sd, method)
  return(10 * log10(sigma.bs))
}

backscatter.xsection.ensemble <- function(scatterer, k, phase.sd=0, n.sim=100) {
  sigma <- rep(0, n.sim)
  for (i in 1:n.sim) {
     sigma[i] <- backscatter.xsection(scatterer, k, phase.sd)
  }
  return(list(mean=mean(sigma), sd=sd(sigma)))
}

#' Backscattering cross section as a function of frequency.
#'
#' @param scatterer A Scatterer object
#' @param freq.start,freq.stop The end points of the frequency range to calculate.
#' @param sound.speed The speed of sound in the surrounding water.
#' @param nfreq The number of frequencies between freq.start and freq.stop.  Defualts
#' to 100.
#' @param ... Additional arguments passed to backscatter.xsection.
#' @return A list with elements freqs, sigma.bs, and TS, containing the frequencies
#' and their corresponding backscatter cross-sections and target strengths.
frequency.spectrum <- function(scatterer, freq.start, freq.stop, sound.speed,
                               nfreq=100, ...) {
  freqs <- seq(freq.start, freq.stop, length.out=nfreq)
  sigma.bs <- rep(0, length(freqs))
  for (i in 1:nfreq) {
    k <- c(0, 0, -2 * pi * freqs[i] / sound.speed)
    sigma.bs[i] <- backscatter.xsection(scatterer, k, ...)
  }
  return(list(freq=freqs, sigma.bs=sigma.bs, TS=10*log10(sigma.bs)))
}

#' Backscattering cross-section as a function of tilt angle.
#'
#' @param scatterer A Scatterer object
#' @param angle.start,angle.stop The end points of the angle range to calculate.
#' @param sound.speed The speed of sound in the surrounding water.
#' @param nangle The number of frequencies between angle.start and angle.stop.  Defualts
#' to 100.
#' @param ... Additional arguments passed to backscatter.xsection.
tilt.spectrum <- function(scatterer, angle.start, angle.stop, freq, sound.speed,
                          nangle=100, ...) {
  angles <- seq(angle.start, angle.stop, length.out=nangle) # degrees
  sigma.bs <- rep(0, length(angles))
  for (i in 1:length(angles)) {
    k <- c(0, 0, -2 * pi * freq / sound.speed)
    sigma.bs[i] <- backscatter.xsection(rotate(scatterer, tilt=angles[i]), k, ...)
  }
  return(list(angle=angles, sigma.bs=sigma.bs, TS=10 * log10(sigma.bs)))
}
