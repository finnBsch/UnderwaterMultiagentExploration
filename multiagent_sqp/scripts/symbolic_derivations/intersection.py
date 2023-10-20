from sympy import *


if __name__ == "__main__":
    x0, y0 = symbols("x0 y0")
    x1, y1 = symbols("x1 y1")
    ox0, oy0 = symbols("ox0 oy0")
    ox1, oy1 = symbols("ox1 oy1")

    k, t = symbols("k t")

    expr0 = x0 + k * (x1 - x0) - (ox0 + t * (ox1 - ox0))
    expr1 = y0 + k * (y1 - y0) - (oy0 + t * (oy1 - oy0))
    res_kt = solve([expr0, expr1], [k, t])
    print(simplify(res_kt[t]))