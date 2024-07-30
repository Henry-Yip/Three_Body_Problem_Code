"""
This code is based on the book Integrable Systems in Celestial Mechanics, Chapter 3. 
We are finding the value of e,  given sigma and R (how to find sigma and R are shown in another code in my repository). 
The intermediate values, c1 and c2, p and a are shown also.  

Author: Henry Yip
"""

import numpy as np

R_list = []
sigma_list = []
time_steps = 100
dt = 0.001
E = -0.5

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
C1 =np.mean(C1_list[:10])
print("C2 Mean is", np.mean(C2_list[:10]))
print("C1 Mean is", np.mean(C1_list[:10]))


if C2 > 0:
    C = np.sqrt(2*C2)
    p = C**2 / 2
    a = -1/E
    e = np.sqrt(1 - p/a)
    print("The value of e is", e)
else:
    C = np.sqrt(2*C1)
    p = C**2 / 2
    a = -1/E
    e = np.sqrt(1 - p/a)
    print("The value of e is", e)
    
