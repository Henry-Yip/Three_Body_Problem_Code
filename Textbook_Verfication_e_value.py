"""
This code is based on the book Integrable Systems in Celestial Mechanics, Chapter 3. 
We are finding p, eta and e for an orbit saved by "python fixed_centres.py --xyz".
The moving body must be the first body in each frame, and the centres must have unit mass and lie at z = 1 and z = -1.
The intermediate values, C1 and C2, are shown also.

Author: Henry Yip
"""

import numpy as np

Input_File_Name = input("Name of input file: ")
dt = float(input("Time between saved points: "))

R_list = []
sigma_list = []

with open(Input_File_Name, 'r') as file:
    lines = file.read().strip().split('\n')

for base_index in range(2, len(lines), 5):
    coords = lines[base_index].split()[1:]
    x, z = float(coords[0]), float(coords[2])
    to_upper = np.hypot(x, z - 1)
    to_lower = np.hypot(x, z + 1)
    R_list.append((to_upper + to_lower)/2)
    sigma_list.append(np.arccos(np.clip((to_lower - to_upper)/2, -1, 1)))

def calculate_c1(R_list, sigma_list, time_steps):
     b = 1
     beta = 0
     mu = 2
     epsilon = 1e-10  # Small value to avoid division by zero

     # Numerical derivative of R with respect to time
     R_dot = np.gradient(R_list, time_steps)
     # Numerical derivative of sigma with respect to time
     sigma_dot = np.gradient(sigma_list, time_steps)
     result = []

     for R, sigma, R_dot_i, sigma_dot_i in zip(R_list, sigma_list, R_dot, sigma_dot):
         cos_sigma = np.cos(sigma)
         cos_sigma_squared = cos_sigma**2

         # Calculate T
         denominator_T1 = R**2 - b**2
         term_T1 = (R**2 - b**2 * cos_sigma_squared) / \
             (denominator_T1 + epsilon) * R_dot_i**2

         term_T2 = (R**2 - b**2 * cos_sigma_squared) * sigma_dot_i**2
         T = 0.5 * (term_T1 + term_T2)

         # Calculate V
         denominator_V = R**2 - b**2 * cos_sigma_squared
         V = -mu * (R + beta * b * cos_sigma) / (denominator_V)

         # Calculate E = T + V
         E = T + V

         # Calculate the terms for the result
         numerator = (R**2 - cos_sigma_squared)**2
         denominator = R**2 - b**2
         term1 = 0.5 * (numerator / (denominator + epsilon)) * R_dot_i**2
         term2 = E * R**2
         term3 = 2 * R

         result.append(term1 - term2 - term3)

     return result


def calculate_c2(R_list, sigma_list, time_steps):
     b = 1
     beta = 0
     mu = 2

     # Numerical derivative of sigma with respect to time
     sigma_dot = np.gradient(sigma_list, time_steps)
     # Numerical derivative of R with respect to time
     R_dot = np.gradient(R_list, time_steps)
     result = []
     term1_list = []
     term2_list = []
     term3_list = []

     for R, sigma, sigma_dot_i, R_dot_i in zip(R_list, sigma_list, sigma_dot, R_dot):
         cos_sigma = np.cos(sigma)
         cos_sigma_squared = cos_sigma**2

         # Calculate T
         T = 0.5 * ((R**2 - b**2 * cos_sigma_squared) / (R**2 - b**2)) * R_dot_i**2 + \
             0.5 * (R**2 - b**2 * cos_sigma_squared) * sigma_dot_i**2

         # Calculate V
         V = -mu * (R + beta * b * cos_sigma) / \
             (R**2 - b**2 * cos_sigma_squared)

         # Calculate E = T + V
         E = T + V

         # Calculate the terms for the result
         term1 = 0.5 * (R**2 - cos_sigma_squared)**2 * sigma_dot_i**2
         term2 = E * cos_sigma_squared
         term3 = E

         term1_list.append(term1)
         term2_list.append(term2)
         term3_list.append(term3)
         result.append(term1 + term2)

     return result, term1_list, term2_list, term3_list  


C2_list, term1_list, term2_list, term3_list = calculate_c2(R_list, sigma_list, dt)
C1_list = calculate_c1(R_list, sigma_list, dt)

C2 = np.mean(C2_list[:10])
C1 = np.mean(C1_list[:10])
E = np.mean(term3_list[:10])  # term3 is the total energy of the moving body
print("C2 Mean is", C2)
print("C1 Mean is", C1)
print("The energy is", E)

# C2 is the constant N in the book and g in the paper.
p = abs(C2)
print("The value of p is", p)
print("The value of eta is", 1/p)

if C2 > 0:
    a = -1/E
    e = np.sqrt(1 - p/a)
    print("The value of a is", a)
    print("The value of e is", e)
else:
    print("C2 is negative (a satellite orbit), so C = sqrt(2*C2) is not real and e is not defined here.")
    
