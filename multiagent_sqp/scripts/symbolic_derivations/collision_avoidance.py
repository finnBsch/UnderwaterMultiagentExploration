
from sympy import *

x0, y0, x1, y1 = symbols("x0, y0, x1, y1")
lagrange_mult = Symbol("lambda")
expr = Matrix([((x1 - x0) ** 2 + (y1 - y0) ** 2) * lagrange_mult])

hess = hessian(expr, [x0, y0, x1, y1, lagrange_mult])
print(hess.eigenvals())  # is PSD
print(Matrix([expr]).jacobian( [x0, y0, x1, y1, lagrange_mult]))
