push!(LOAD_PATH, expanduser("~/Documents/nano-optics/CoupledDipole.jl/"))
using Revise
using CoupledDipole

# depolarisation <- function (x1, x2 = x1, x3 = x2)
# {
#   ## scaled version to help integration
#   V = x1*x2*x3

#   integrand <- function(q, r, s) {
#     1/((1 + q) * sqrt((q + 1) * (q + r^2) * (q + s^2)))
#   }

#   I1 <- integrate(integrand, r=x2/x1, s=x3/x1, lower = 0, upper = Inf)
#   I2 <- integrate(integrand, r=x1/x2, s=x3/x2, lower = 0, upper = Inf)
#   I3 <- integrate(integrand, r=x1/x3, s=x2/x3, lower = 0, upper = Inf)

#   V/2 * c(I1$value / x1^3,
#           I2$value / x2^3,
#           I3$value / x3^3)
# }


function integrand(q, r, s)
  1.0 / ((1.0 + q) * sqrt((q + 1.0) * (q + r^2) * (q + s^2)))
end

using QuadGK
a = quadgk(q -> integrand(q, 0.1, 0.1), 0.0, Inf)

depolarisation_ellipsoid(40.01, 40.0, 40.02)


depolarisation_ellipsoid(20.01, 20.0, 50.02)
depolarisation_ellipsoid(20.0, 20.01, 50.0)
depolarisation_ellipsoid(50.01, 20.0, 20.0)
depolarisation_ellipsoid(20.01, 50.0, 20.0)
