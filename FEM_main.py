import FEM
import numpy as np
import main_programv2

gridsize = 0.125
hstep = gridsize
threshold = 0.01
lbda = 0.1
mu = 3
m = FEM.mesh(gridsize, threshold)
m.plot()
s = FEM.stiffness_matrix(m, lbda, mu)
q = np.zeros((2 * len(m.points), 1))
farr = [main_programv2.linear]
methodname = 'Shortest path fit with Euler forward relaxation'
newtriang = False
flip = False
simplepointrelax = False
simplegridtozeropoint = False
fixedpointrelaxation = False
eulerrelaxation = False
shortestpath = True
distplot = True
redistribute = True
m.initzeroinfo(farr[0])
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

u = FEM.meshrelaxation(m, u_xy, lbda, mu)

for p in m.points:
    p.coordinates[0] += float(u[p.index])
    p.coordinates[1] += float(u[p.index + len(m.points)])
m.plot()
m.setupboundaryplist()
m.finallevelsetinfoupdate()
for i in m.points:
    m.points[i.index].iszeropath = m_fit.points[i.index].iszeropath
for i in range(len(m.edges)):
    m.edges[i].islevelset = m_fit.edges[i].islevelset
for e in m.edges:  # without altering structure (no level-set lines are deleted)
    e.length = e.calclength()
for el in m.elements:
    el.size = el.calcsize()
    el.skewness = el.calcskew()
m.plotzeroz('testfigures/test.png', mode=2)
m_fit.plotzeroz('testfigures/test2.png', mode=2)
