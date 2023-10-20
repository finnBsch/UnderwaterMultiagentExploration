import numpy as np

file = open("hess.txt", "r")
lines = file.readlines()
file.close()

header = lines[0]
print(header.split())
mat_size = header.split()[4]
temp_ = mat_size.split("-")

nx = int(temp_[0])

ny = int(''.join(c for c in temp_[2] if c.isdigit()))
mat = np.zeros((nx, ny), dtype=np.float32)
lines = lines[1:]
for line in lines:
    els = line.split()
    id_x = int(''.join(c for c in els[0] if c.isdigit()))
    id_y = int(''.join(c for c in els[1] if c.isdigit()))
    val = float(els[3])
    mat[id_x, id_y] = val


eigs =np.linalg.eig(mat)[0]
print(eigs)
print(np.all(eigs > -1e-3))
