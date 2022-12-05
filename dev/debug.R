library(cda)
# library(rgl)
library(ggplot2)
library(dielectric)
library(cubs)
library(tidyr)
library(dplyr)
library(purrr)
library(progress)


wavelength = c(400, 500)
gold = epsAu(wavelength)

cl = structure(list(positions = structure(c(-50, 0, 0, 50, 0, 00), dim = c(3,2)),
                    sizes = structure(c(50, 20, 20, 50, 20, 20), dim = c(3,2)),
                    angles = structure(c(0, 0, 0, 0, 0, 0), dim = c(3,2))),
               class = "cluster")

cl0 = structure(list(positions = structure(c(0, 0, 0), dim = c(3,1)),
                    sizes = structure(c(50, 20, 20), dim = c(3,1)),
                    angles = structure(c(0, 0, 0), dim = c(3,1))),
               class = "cluster")

cl0
cl$positions


ii = 1;

lambda = gold$wavelength[ii]
n_medium = 1.33
kn = n_medium * 2*pi / lambda

# gold[ii,]
# -1.6496568840720354 + 5.7717630808981655im

# alpha_kuwata(500, -10+1im, 1.33^2, SVector(30, 30, 50))
# 3-element SVector{3, ComplexF64} with indices SOneTo(3):
#   77076.04648078185 + 26235.664281642246im
# 77076.04648078185 + 26235.664281642246im
# -98187.15974124733 + 205835.30299929058im
# V <- 4 * pi/3 * 30*30*50
# chi <- cda:::depolarisation(30,30,50)
# cda:::alpha_kuwata(500, -10+1i, V, 30, chi[1], 1.33)
# cda:::alpha_kuwata(500, -10+1i, V, 50, chi[3], 1.33)

Alpha = cda::alpha_ellipsoid(sizes = cl0$sizes, material = gold[ii,], medium = n_medium)

# Alpha = alpha_particles(λ, ε, n_medium^2, cl.sizes)
# 2-element Vector{SVector{3, ComplexF64}}:
# [10898.262981422746 + 18516.11377073886im, -903.6988080340118 + 24897.841185369187im, 10898.262981422744 + 18516.113770738863im]
# [10898.262981422746 + 18516.11377073886im, -903.6988080340118 + 24897.841185369187im, 10898.262981422744 + 18516.113770738863im]

# alpha_kuwata(400, ε, 1.33^2, SVector(20, 20, 50))
# 3-element SVector{3, ComplexF64} with indices SOneTo(3):
#   12403.616901231882 + 13704.077043611223im
# 12403.616901231882 + 13704.077043611223im
# -6869.361725855006 + 22867.056384367785im

ParticleRotations = cl.rotations
AlphaBlocks = map((R, A) -> R' * diagm(A) * R, ParticleRotations, Alpha)



