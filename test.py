from main_programv2 import *

m = mesh(hstep, hstep * 0.3)
m.plot()

hstep = 0.125
farr = [circle]
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
m_fit = simplemeshfit(0.3 * hstep, methodname, hstep, hstep / 100, farr, newtriang, flip, distplot, simplepointrelax,
                      simplegridtozeropoint, redistribute)
m_fit.plot()
u_x_fit = []
u_y_fit = []
for i in range(len(m.points)):
    u_x_fit.append(m.points[i].coordinates[0] - m_fit.points[i].coordinates[0])
    u_y_fit.append(m.points[i].coordinates[1] - m_fit.points[i].coordinates[1])
