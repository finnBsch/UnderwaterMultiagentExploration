import sympy as sy
import numpy as np

xk, yk, phik, vk = sy.symbols("xk, yk, phik, vk")
xk1, yk1, phik1, vk1 = sy.symbols("xk1, yk1, phik1, vk1")

dxk, dyk, dphik, dvk = sy.symbols("dxk, dyk, dphik, dvk")
dxk1, dyk1, dphik1, dvk1 = sy.symbols("dxk1, dyk1, dphik1, dvk1")

u1, u2 = sy.symbols("u0, u1")
du1, du2 = sy.symbols("du0, du1")

dt = sy.Symbol("dt")

xmuk, ymuk, phimuk, vmuk = sy.symbols("xmuk, ymuk, phimuk, vmuk")
dxmuk, dymuk, dphimuk, dvmuk = sy.symbols("dxmuk, dymuk, dphimuk, dvmuk")

muk_v = [xmuk, ymuk, phimuk, vmuk]
g1 = xk1 - xk - dt*sy.cos(phik)*vk
g2 = yk1 - yk - dt*sy.sin(phik)*vk
g3 = phik1 - phik - dt*u1
g4 = vk1 - vk - dt*u2

g = sy.Matrix([g1, g2, g3, g4])

mult_vec = sy.Matrix([dxk, dyk, dphik, dvk, dxk1, dyk1, dphik1, dvk1, du1, du2])

dg = g.jacobian([xk, yk, phik, vk, xk1, yk1, phik1, vk1, u1, u2])

print("dg * dz")
print(dg@mult_vec)

mult_vec = sy.Matrix([dxk, dyk, dphik, dvk, dxk1, dyk1, dphik1, dvk1, du1, du2, dxmuk, dymuk, dphimuk, dvmuk])
for i in range(4):
    gg = sy.Matrix([muk_v[i] * g[i]])
    print(i)
    ddg = sy.hessian(gg, [xk, yk, phik, vk, xk1, yk1, phik1, vk1, u1, u2, xmuk, ymuk, phimuk, vmuk])
    # print(sy.simplify(ddg))
    print(sy.simplify(mult_vec.T @ ddg @ mult_vec))
