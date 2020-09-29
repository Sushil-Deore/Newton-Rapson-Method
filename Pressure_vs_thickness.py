# importing important modules

import matplotlib.pyplot as plt
import numpy as np
import math


# Defining function of Mullar's Equation

def func(p, sigma, beta, r, h):
    # spliting equcation into terms
    term_1 = pow(p, 3) * (1 - pow(beta, 2))

    term_2 = (0.4 * h * pow(beta, 2) - (sigma * pow(h, 2) / pow(r, 2))) * pow(p, 2)

    term_3 = (pow(sigma, 2) * pow(h, 4) / (3 * pow(r, 4))) * p

    term_4 = pow((sigma * pow(h, 2) / (3 * pow(r, 2))), 3)

    return term_1 + term_2 + term_3 - term_4


# Defining function of Mullar's Equation after differentiating

def func_diff(p, sigma, beta, r, h):
    # spliting equcation into terms
    diff_t1 = 3 * pow(p, 2) * (1 - pow(beta, 2))

    diff_t2 = (0.4 * h * pow(beta, 2) - (sigma * pow(h, 2) / pow(r, 2))) * 2 * p

    diff_t3 = (pow(sigma, 2) * pow(h, 4) / (3 * pow(r, 4)))

    return diff_t1 + diff_t2 + diff_t3


# Input Parameters

alpha = 1
sigma = 150  # pounds per square inch (psi)
beta = 0.5
r = 40  # feet
h = [0.6, 1.2, 1.8, 2.4, 3, 3.6, 4.2]  # feet
tol = 1e-14
iter = 1

# Performing newton iteration

p_guess = 120  # to change
pressure = []

# Defining Print function for Title of Thickness & Pressure
print("Thickness  Pressure")

for i in range(0, len(h)):
    height = h[i]

    while (abs(func(p_guess, sigma, beta, r, height)) > tol):
        p_guess = p_guess - alpha * (func(p_guess, sigma, beta, r, height) / func_diff(p_guess, sigma, beta, r, height))
        iter = iter + 1

    pressure.append(p_guess)

    print('%f   %f' % (height, pressure[i]))
