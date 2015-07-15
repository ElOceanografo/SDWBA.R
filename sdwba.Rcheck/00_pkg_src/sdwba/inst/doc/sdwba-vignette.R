## ------------------------------------------------------------------------
library(sdwba)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(generic.krill.Conti2006, format="markdown")

## ------------------------------------------------------------------------
generic.krill.m <- generic.krill.Conti2006
generic.krill.m[c("x", "y", "z", "a")] <- generic.krill.m[c("x", "y", "z", "a")] * 1e-3
krill <- with(generic.krill.m, Scatterer(x, y, z, a, g, h))

## ---- eval=FALSE---------------------------------------------------------
#  scat <- load.scatterer("path/to/my_scatterer.csv")

## ---- fig.width=7--------------------------------------------------------
plot(krill)

## ------------------------------------------------------------------------
sound.speed.water <- 1480
# density.water <- 1027
k.mag <- 120e3 * 2 * pi / sound.speed.water
k <- c(0, 0, -k.mag)

## ---- fig.width=7, fig.height=7------------------------------------------
krill.2 <- rotate(krill, tilt=45)
plot(krill.2)

## ---- fig.width=7, fig.height=7------------------------------------------
backscatter.xsection(krill, k)
target.strength(krill, k)

## ---- fig.width=7, fig.height=5------------------------------------------
n.rep <- 100
sigma.bs <- rep(0, n.rep)

for (i in 1:n.rep) {
  sigma.bs[i] <- backscatter.xsection(krill.2, k, phase.sd=10)
}
hist(10 * log10(sigma.bs), 20, main="TS histogram")

10 * log10(mean(sigma.bs))
sd(10 * log10(sigma.bs))

## ---- fig.width=7, fig.height=5------------------------------------------
freq.spec <- frequency.spectrum(krill, 12e3, 400e3, sound.speed.water, nfreq=120)
plot(freq.spec$freq, freq.spec$TS, ty='l')

tilt.spec <- tilt.spectrum(krill, -180, 180, freq=120e3, sound.speed=sound.speed.water, nangle=180)
plot(tilt.spec$angle, tilt.spec$TS, ty='l')

## ------------------------------------------------------------------------
# fs.cont <- frequency.spectrum(krill, 12e3, 400e3, sound.speed.water,
#                          nfreq=120, method="continuous")
# plot(TS ~ freq, fs, ty='l')
# lines(TS ~ freq, fs.cont, col='blue')

