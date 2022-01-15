# Palatini_TOV

This repo contains the code for my master thesis with the title 'Palatini f(R,Q)-gravity applied on the structure of neutron stars'
You can find the pdf of the thesis here: https://unipub.uni-graz.at/obvugrhs/download/pdf/2679877?originalFilename=true 

## Background
Stellar structure equations allow to calculate the mass of spherical symmetric objects (such as neutron stars), by integrating over interdependent pressure and mass gradients.

The gradients and therefor the total mass of a star crucially depend on two factors
- the Equation of State (EOS), characterizing the relation between density and pressure within a star
- the underlying theory of gravity

By chosing a set of EOS and theory of gravity one can compute the total mass of neutron stars, given the density in the center of the star.
The computed total mass profile for different center densities can then be compared to the observed masses of neutron stars and give insights
into the underlying EOS of stars and the validity of theories of gravity.

## What the code does
The code implements the numerical integration over density of a spherical symmetric object via the 4th-order Runge-Kutta method (RK4).
As parameters it takes the EOS either in tabular form or as analytical expression as well as the stellar structure equation from a theory of gravity.
In this version I looked into the results using 1) Einsteins theory of gravity and 2) modified theories of gravity (f(R), f(R,Q))
f(R) and f(R,Q) theories extend the action S in the Einstein formulation with additional scalars built from the Ricci Scalar R or the quadratic curvature invariant Q.

![Screenshot from 2022-01-15 17-21-17](https://user-images.githubusercontent.com/38614357/149629321-e0b83df1-9ff0-443b-ac54-7a0aa3903682.png)

