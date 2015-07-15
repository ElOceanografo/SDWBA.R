context("Accuracy")

sound.speed <- 1456

generic.krill.m <- generic.krill.Conti2006
generic.krill.m[c("x", "y", "z", "a")] <- generic.krill.m[c("x", "y", "z", "a")] * 1e-3
krill <- with(generic.krill.m, Scatterer(x, y, z, a, g, h))


freqs <- c(38, 120, 200, 333) * 1e3

for (f in freqs) {
  k <- c(0, 0, -f) * 2 * pi / sound.speed
  print(target.strength(krill, k))
}
