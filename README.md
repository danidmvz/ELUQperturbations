# On the propagation of uncertainty in Eulerian-Lagrangian point-particle methods with non-deterministic forcing laws

An analysis of the propagation of uncertainty in the particle phase solution of particle-laden flows described with the
point-particle approach is presented by combining a regular perturbation problem with the method of moments. Non-
deterministic forcing laws account for uncertainties in the point-particles’ force related to the modeling of the near
flow, particle-particle collisions, particle shapes and other effects that are inaccessible by general deterministic forcing
laws with a reduced number of parameters. These uncertainties are recasted in a random forcing model with random
coefficients that do not depend on time. A regular perturbation problem can be formulated to solve the particle equations
for the case in which such random coefficients present small fluctuations around the average value. Combining this
with the method of moments, explicit expressions for the moments of the particle phase in closed form are found. Also,
Monte Carlo simulations can be highly accelerated. Samples can be computed by simply a linear combination of the
successive orders of the perturbation problem that are solved only once, instead of solving the original point-particle
equations for each sample. Examples of this approach are presented and validated with Monte Carlo simulations for a
variety of test cases.