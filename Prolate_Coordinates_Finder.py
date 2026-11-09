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

# R is the mean of the distances to the two centres, at z = 1 and z = -1.
# Unlike R = z/cos(sigma), this stays accurate when z is very small.
R = (math.hypot(x, z - 1) + math.hypot(x, z + 1))/2

print(f"The value of R is: {R}")
print(f"The value of σ is: {sigma}")