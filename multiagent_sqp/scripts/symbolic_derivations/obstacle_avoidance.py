from sympy import *
from sympy.printing import cxxcode

obs_x0, obs_y0, obs_x1, obs_y1 = symbols("obs_x0, obs_y0, obs_x1, obs_y1")

x, y = symbols("x, y")

dx, dy = symbols("dx, dy")

length_obs = (obs_x1 - obs_x0)**2 + (obs_y1 - obs_y0)**2
t = ((x - obs_x0) * (obs_x1 - obs_x0) + (y - obs_y0) * (obs_y1 - obs_y0))/length_obs
px = obs_x0 + t * (obs_x1 - obs_x0)
py = obs_y0 + t * (obs_y1 - obs_y0)
distance = (px - x)**2 + (py - y)**2

expr = Matrix([distance])
jac = expr.jacobian([x, y])
print("x : ")
print(simplify(jac[0, 0]))
print("y : ")
print(simplify(jac[0, 1]))

print("all : ")
print(cxxcode(simplify(jac@Matrix([dx, dy]))[0, 0]))
