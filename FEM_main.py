import FEM
import numpy as np
from main_programv2 import simplemeshfit, dumbbell, linear, circle, star, cos_times_sin

gridsize = 0.125
hstep = gridsize
threshold = 0.01
lbda = 0.05
mu = 1
m = FEM.mesh(gridsize, threshold)
m.plot()
s = FEM.stiffness_matrix(m, lbda, mu)
q = np.zeros((2 * len(m.points), 1))
farr = [linear]
methodname = 'FEM'
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
m_fit = simplemeshfit(0.3 * hstep, methodname, hstep, hstep / 100, farr, newtriang, flip, distplot,
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
    m.edges[i].iszeroedge = m_fit.edges[i].iszeroedge
for e in m.edges:
    e.length = e.calclength()
for el in m.elements:
    el.size = el.calcsize()
    el.skewness = el.calcskew()

m.plotzeroz('figures_FEM/FEMfit.png', mode=2)
m_fit.plotzeroz('figures_FEM/FEMstartsituation.png', mode=2)

redistribute = False
eulerrelaxation = True
m_efr = simplemeshfit(0.3 * hstep, 'Shortest path fit with Euler forward relaxation', hstep, hstep / 100, farr,
                      newtriang, flip, distplot, simplepointrelax, simplegridtozeropoint, redistribute)
femstats, efrstats = FEM.qualitycomparison(m, m_efr)
FEM.qualitycomparisonplot(femstats, efrstats, lbda, mu, farr[0].__name__)
