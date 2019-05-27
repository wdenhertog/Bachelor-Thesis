import FEM
import numpy as np
import main_programv2

gridsize = 0.125
hstep = 0.125
threshold = 0.01
lbda = 5
mu = 3
m = FEM.mesh(gridsize, threshold)
m.plot()
s = FEM.stiffness_matrix(m, lbda, mu)
q = np.zeros((2 * len(m.points), 1))
farr = [main_programv2.star]
methodname = 'Shortest path fit with redistribution'
newtriang = False
flip = False
simplepointrelax = False
simplegridtozeropoint = False
fixedpointrelaxation = False
eulerrelaxation = False
shortestpath = True
distplot = True
redistribute = True
m_fit = main_programv2.simplemeshfit(0.3 * hstep, methodname, hstep, hstep / 100, farr, newtriang, flip, distplot,
                                     simplepointrelax,
                                     simplegridtozeropoint, redistribute)
m_fit.plot()
u_x_fit = np.zeros(len(m.points))
u_y_fit = np.zeros(len(m.points))
for i in range(len(m.points)):
    u_x_fit[i] = (m_fit.points[i].coordinates[0] - m.points[i].coordinates[0])
    u_y_fit[i] = (m_fit.points[i].coordinates[1] - m.points[i].coordinates[1])
u_xy = np.hstack((u_x_fit, u_y_fit))
# for p in m.boundarypoints:
#     q[378] = -0.5
#     q[377] = -0.5
#     q[376] = -0.9
#     if p.coordinates[0] == -1 or p.coordinates[0] == 1:
#         for i in range(2*len(m.points)):
#             s[p.index][i] = 0
#         s[p.index][p.index] = 1
#     if p.coordinates[1] == -1:
#         for i in range(2 * len(m.points)):
#             s[p.index][i] = 0
#             s[p.index + len(m.points)][i] = 0
#         s[p.index][p.index] = 1
#         s[p.index + len(m.points)][p.index + len(m.points)] = 1

u = FEM.meshrelaxation(m, u_xy, lbda, mu)

for p in m.points:
    p.coordinates[0] += float(u[p.index])
    p.coordinates[1] += float(u[p.index + len(m.points)])
m.plot()
m.plotzeroz('testfigures/test.png', mode=2)
