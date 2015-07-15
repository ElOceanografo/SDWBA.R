
context("Function execution")

sound.speed.water <- 1480
density.water <- 1027

generic.krill.m <- generic.krill.Conti2006
generic.krill.m[c("x", "y", "z", "a")] <- generic.krill.m[c("x", "y", "z", "a")] * 1e-3
krill <- with(generic.krill.m, Scatterer(x, y, z, a, g, h))

save.scatterer(krill, filename = "test_scatterer.csv")
krill2 <- load.scatterer("test_scatterer.csv")
file.remove("test_scatterer.csv")

plot(krill)
file.remove("Rplots.pdf")

rotate(krill, 45, 45, 45)
rescale(krill, 2, 2, 2, 2, 2)

k.mag <- 120e3 * 2 * pi / sound.speed.water
k <- c(0, 0, -k.mag)
backscatter.xsection(krill, k, phase.sd=10)
target.strength(krill, k, phase.sd=10)
k <- c(0, 0, -k.mag)
backscatter.xsection(krill, k, phase.sd=10, method="continuous")
target.strength(krill, k, phase.sd=10, method="continuous")

freq.spec <- frequency.spectrum(krill, 12e3, 400e3, sound.speed.water, nfreq=12)

tilt.spec <- tilt.spectrum(krill, -180, 180, freq=120e3, sound.speed=sound.speed.water, nangle=18)
