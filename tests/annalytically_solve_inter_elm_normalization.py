"""Analytically normalize the inter-ELM function of the form:


    q_parallel = {R exp(-(R - Rbdry)/ lambda_n), for Rbdry + Rbreak > R > Rbdry
                 {R exp(-(R - Rbdry) / lambda_m + Rbreak / lambda_m - Rbreak / lambda_n), for Rwall > R > Rbdry + Rbreak

With this we imply that the function does not spread towards infinity but
infact stops at an end.

"""
from sympy import exp, Symbol, integrate
from sympy.physics.units.quantities import Quantity


def pprint(sol):
    x = str(sol)
    x = x.replace("LAMBDA_NEAR", "ln")
    x = x.replace("LAMBDA_MAIN", "lm")
    x = x.replace("R_BOUNDARY", "R_bdry")
    x = x.replace("R_WALL", "R_wall")
    x = x.replace("R_BREAK", "R_b")
    x = x.replace("exp", "np.exp")
    print(x)
R = Symbol('R')
ln = Quantity("LAMBDA_NEAR")
lm = Quantity("LAMBDA_MAIN")

Rbdry = Quantity("R_BOUNDARY")
Rwall = Quantity("R_WALL")
# Rwall = Infinity()
Rbreak = Quantity("R_BREAK")

A = Rbdry
B = Rbdry + Rbreak
C = Rwall

# From Rbdry to Rbdry + Rbreak
sol_near = integrate(R * exp(-(R-Rbdry)/ln), (R, A, B))

# From Rbdry + Rbreak to Rwall
sol_main = integrate(R * exp((-(R-Rbdry)/lm) + (Rbreak/lm) - (Rbreak/ln)), (R, B, C))

# pprint(sol_near)
# pprint(sol_main)

sol = sol_near + sol_main

pprint(sol_near.simplify())
pprint(sol_main.simplify())
final_sol = sol
pprint(final_sol)

