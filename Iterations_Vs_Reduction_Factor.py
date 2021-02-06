# importing important modules

import matplotlib.pyplot as plt
import numpy as np
import math
import warnings
warnings.filterwarnings('ignore')


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

sigma = 150  # pounds per square inch (psi)
beta = 0.5
r = 40  # feet
h = 0.6  # feet
tol = 1e-2

# Performing newton iteration

alpha = []
i = 0.1
j = 0

while (i < 5.5):
    alpha.append(round(i, 2))
    i = i + 0.1
    j = j + 1

pressure = []
iter = []

for i in range(0, len(alpha)):
    p_guess = 120
    iter_1 = 1

    while (abs(func(p_guess, sigma, beta, r, h)) > tol):
        p_guess = p_guess - alpha[i] * (func(p_guess, sigma, beta, r, h) / func_diff(p_guess, sigma, beta, r, h))

        iter_1 = iter_1 + 1

    iter.append(iter_1 - 1)
    pressure.append(p_guess)

min_iter = min(iter)

for i in range(0, len(iter)):
    if (min_iter == iter[i]):
        pos = i
        break

print(f"The Optimum RF is {alpha[i]}")

# Plotting No of iterations vs. Reduction Factor

plt.plot(alpha, iter)

plt.xlabel('No. of Iterations')
plt.ylabel('Reduction Factor')

plt.title('No. of Iterations Vs Reduction Factor')
plt.plot(alpha[i], min_iter, 'o')

plt.legend(['Optimum RF'])
plt.show()

