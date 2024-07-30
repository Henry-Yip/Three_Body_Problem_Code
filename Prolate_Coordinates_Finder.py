"""
The code calculates the planar prolate spheroidal coordinates given the necessary cartesian coordinates.
The derivation is given here: https://henry-yip.github.io/EulerProblemPlots/
This code is based on the description in "Integrable Systems in Celestial Mechanics" 

Author: Henry Yip
"""
import math
import numpy as np

x = float(input("Please enter the x-coordinate: "))
z = float(input("Please enter the z-coordinate: "))

def calculate_sigma(x, z):
    x = abs(x)  # Use the absolute value of x for the calculation
    u = (-(x**2 + z**2 - 1) + np.sqrt((x**2 + z**2 - 1)**2 + 4 * x**2)) / 2
    u_sqrt = np.sqrt(u)
    sigma = np.arcsin(u_sqrt)
    if z < 0:
        return np.pi - sigma

    else:
        return sigma
    
sigma = calculate_sigma(x,z)

if z != 0:
     R = z/(math.cos(sigma))
else:
     sin_sigma = math.sin(sigma)

     R_squared = (x/sin_sigma)**2+1

     R = math.sqrt(R_squared)

print(f"The value of R is: {R}")
print(f"The value of σ is: {sigma}")