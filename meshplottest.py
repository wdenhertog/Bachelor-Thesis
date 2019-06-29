import FEM
import numpy as np
import main_programv2

gridsize = 0.125
hstep = 0.125
threshold = 0.01
lbda = 5
mu = 3
m = main_programv2.mesh(gridsize, threshold)
m.plot()
s = FEM.stiffness_matrix(m, lbda, mu)
q = np.zeros((2 * len(m.points), 1))
farr = [main_programv2.circle]
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
m.plotzeroz('testfigures/test.png', mode=2)