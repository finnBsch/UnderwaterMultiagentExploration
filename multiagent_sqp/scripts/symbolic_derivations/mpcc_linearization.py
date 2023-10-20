"""
MPCC local linearization for arbitrary reference path parameterized by progress variable theta.
"""

import numpy as np
from sympy import * 

x, y, theta = symbols("x, y, theta")
alpha = Function("alpha")(theta)
x_ = Function("x_")(theta)
y_ = Function("y_")(theta)

e_c = sin(alpha) * (x - x_) - cos(alpha) * (y - y_)
e_l = -cos(alpha) * (x - x_) - sin(alpha) * (y - y_)

dec_sq = Matrix([e_c**2]).jacobian([x, y, theta])
del_sq = Matrix([e_l**2]).jacobian([x, y, theta])
print("Contouring cost gradients: ")
print("d/dx:")
print(simplify(dec_sq[0]))
print("d/dy:")
print(simplify(dec_sq[1]))
print("d/dtheta:")
print(simplify(dec_sq[2]))
print("")
print("Contouring cost gradients: ")
print("d/dx:")
print(simplify(del_sq[0]))
print("d/dy:")
print(simplify(del_sq[1]))
print("d/dtheta:")
print(simplify(del_sq[2]))
