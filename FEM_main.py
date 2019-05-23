import FEM
import numpy as np

gridsize = 1
threshold = 0.01
lbda = 5
mu = 3
m = FEM.mesh(gridsize, threshold)
m.plot()
s = FEM.stiffness_matrix(m, lbda, mu)
q = np.zeros((2 * len(m.points), 1))
for p in m.boundarypoints:
    if p.coordinates[1] == 1:
        q[p.index + len(m.points)] = -0.5
    elif p.coordinates[1] == -1:
        for i in range(2 * len(m.points)):
            s[p.index][i] = 0
            s[p.index + len(m.points)][i] = 0
        s[p.index][p.index] = 1
        s[p.index + len(m.points)][p.index + len(m.points)] = 1

u = np.linalg.solve(s, q)

for p in m.points:
    p.coordinates[0] += float(u[p.index])
    p.coordinates[1] += float(u[p.index + len(m.points)])
m.plot()
