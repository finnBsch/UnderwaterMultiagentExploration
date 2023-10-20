"""
CENTER OF MASS APPROACH
MPCC local linearization for piecewise linear track (constructed from linear segments). Each track segment starts at
 x0, y0 and has angle alpha. Thus, we can write (y - y0)/(x - x0) = tan(alpha). Consequently, we can remove the DOF y
 by writing y = sin(alpha)/cos(alpha) * (x - x0) + y0. In the script, this means that y_ can be replaced by
 sin(alpha)/cos(alpha) * (x_ - x0) + y0
"""

import numpy as np
from sympy import *

x_a, y_a, theta = symbols("x_a, y_a, theta")
x_b, y_b = symbols("x_b, y_b")

x = (x_a + x_b)/2
y = (y_a + y_b)/2
dx_a, dy_a, dx_b, dy_b, dtheta = symbols("dx_a, dy_a, dx_b, dy_b, dtheta")
x0, y0 = symbols("x0, y0")
alpha = Symbol("alpha")  # Constant over theta
x_ = Function("x_")(theta)
# y_ = Function("y_")(theta) replaced!
y_ = tan(alpha) * (x_ - x0) + y0

dvars = Matrix([dx_a, dy_a, dx_b, dy_b, dtheta])
e_c = simplify(sin(alpha) * (x - x_) - cos(alpha) * (y - y_))
e_l = simplify(-cos(alpha) * (x - x_) - sin(alpha) * (y - y_))
print("Costs: ")
print("Contouring squared: ")
print(simplify(pow(e_c, 2)))
print("Lag squared: ")
print(simplify(pow(e_l, 2)))
dec_sq = Matrix([e_c**2]).jacobian([x_a, y_a, x_b, y_b, theta])
del_sq = Matrix([e_l**2]).jacobian([x_a, y_a, x_b, y_b, theta])
print(dec_sq)
print("Contouring cost gradients: ")
print("d/dx_a:")
print(simplify(dec_sq[0]))
print("d/dy_a:")
print(simplify(dec_sq[1]))
print("d/dx_b:")
print(simplify(dec_sq[2]))
print("d/dy_b:")
print(simplify(dec_sq[3]))
print("d/dtheta:")
print(simplify(dec_sq[4]))

print("")
print("Lag cost gradients: ")
print("d/dx_a:")
print(simplify(del_sq[0]))
print("d/dy_a:")
print(simplify(del_sq[1]))
print("d/dx_b:")
print(simplify(del_sq[2]))
print("d/dy_b:")
print(simplify(del_sq[3]))
print("d/dtheta:")
print(simplify(del_sq[4]))

Hess_c = hessian(e_c**2, (x_a, y_a, x_b, y_b, theta))
print("Quadratic terms: ")
print("Contouring: ")
print(simplify(dvars.transpose() @ Hess_c @ dvars))
Hess_l = hessian(e_l**2, (x_a, y_a, x_b, y_b, theta))
print("Lag: ")
t = simplify(dvars.transpose() @ Hess_l @ dvars)
t = simplify(t.subs(Derivative(x_, theta), cos(alpha)))
print(t)


