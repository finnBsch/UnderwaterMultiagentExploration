"""
MPCC local linearization for piecewise linear track (constructed from linear segments). Each track segment starts at
 x0, y0 and has angle alpha. Thus, we can write (y - y0)/(x - x0) = tan(alpha). Consequently, we can remove the DOF y
 by writing y = sin(alpha)/cos(alpha) * (x - x0) + y0. In the script, this means that y_ can be replaced by
 sin(alpha)/cos(alpha) * (x_ - x0) + y0
"""

import numpy as np
from sympy import *

x, y, theta = symbols("x, y, theta")
dx, dy, dtheta = symbols("dx, dy, dtheta")
x0, y0 = symbols("x0, y0")
alpha = Symbol("alpha")  # Constant over theta
x_ = Function("x_")(theta)
# y_ = Function("y_")(theta) replaced!
y_ = tan(alpha) * (x_ - x0) + y0

dvars = Matrix([dx, dy, dtheta])
e_c = simplify(sin(alpha) * (x - x_) - cos(alpha) * (y - y_))
e_l = simplify(-cos(alpha) * (x - x_) - sin(alpha) * (y - y_))
print("Error Lag: ")
print(simplify(e_l))
print("Error Contouring: ")
print(simplify(e_c))
print("Costs: ")
print("Contouring squared: ")
print(simplify(pow(e_c, 2)))
print("Lag squared: ")
print(simplify(pow(e_l, 2)))
dec_sq = Matrix([e_c**2]).jacobian([x, y, theta])
del_sq = Matrix([e_l**2]).jacobian([x, y, theta])
print("Contouring cost gradients: ")
print("d/dx:")
print(simplify(dec_sq[0])) # ((x_a + x_b - 2*x_(theta))*sin(alpha) - (-2*y0 + y_a + y_b + 2*(x0 - x_(theta))*tan(alpha))*cos(alpha))*sin(alpha)/2
print("d/dy:")
print(simplify(dec_sq[1]))
print("d/dtheta:")
print(simplify(dec_sq[2]))
print("")
print("Lag cost gradients: ")
print("d/dx:")
print(simplify(del_sq[0]))
print("d/dy:")
print(simplify(del_sq[1]))
print("d/dtheta:")
print(simplify(del_sq[2]))

Hess_c = hessian(e_c**2, (x, y, theta))
print("Quadratic terms: ")
print("Contouring: ")
print(simplify(dvars.transpose() @ Hess_c @ dvars))
Hess_l = hessian(e_l**2, (x, y, theta))
print("Lag: ")
t = simplify(dvars.transpose() @ Hess_l @ dvars)
t = simplify(t.subs(Derivative(x_, theta), cos(alpha)))
print(t)


