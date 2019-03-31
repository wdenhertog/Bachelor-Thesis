# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 13:04:47 2018

@author: douden
"""

#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as pyplot
import matplotlib.tri as tri
import matlab
import matlab.engine as me
from mpl_toolkits.mplot3d import Axes3D



def create_mesh(x,y,h):
    matlab_engine= me.start_matlab()
    x = matlab_engine.eval(np.array2string(x))
    y = matlab_engine.eval(np.array2string(y))
    h = matlab_engine.eval(np.array2string(h))
    matlab_engine.workspace['x'] = x
    matlab_engine.workspace['y'] = y
    matlab_engine.workspace['h'] = h
    matlab_engine.mesh_gen(nargout=0)
    elements = matlab_engine.workspace['elements']
    elements = np.array(elements).astype(int)
    edges = matlab_engine.workspace['edges']
    edges = np.array(edges).astype(int)
    points = matlab_engine.workspace['points']
    points = np.array(points)        
    matlab_engine.quit()
    return elements, edges, points


def findedgesindices(edges, edgelist):
    edges.sort()
    indices = []
    for e in edges:
        i = e[0]
        while edgelist[i][0] < e[0]:
            i += e[0] - edgelist[i][0]
        while edgelist[i][0] == e[0]:
            if edgelist[i][1] == e[1]:
                break
            i += 1
        indices.append(i)
    return indices
    


class meshclass:
    def __init__(self, hstep,
                     x = np.array([-1,1,1,-1]),
                     y = np.array([-1,-1,1,1])   ):
        h = np.array(hstep)
            # check whether saved mesh exists
        file_str = "h"+np.array2string(h)+"_x"+np.array2string(x)+"_y"+np.array2string(y)+".npz"
        try:
            npzfile = np.load(file_str)
            print("Mesh does exist\nLoading the mesh")
            for key in npzfile.keys():
                exec("self."+ key + " = npzfile[key]")
            print("Mesh loaded")
        except IOError:
            print("Mesh does not exist\nCreating the mesh")
            # create mesh, random numbers
            self.elements, self.edges, self.points = create_mesh(x,y,h)
            print("Saving the mesh")
            np.savez(file_str, elements=elements,edges=edges,points=points)
            print("Mesh saved")

    def plot(self):
        triang = tri.Triangulation(self.points[0,:], self.points[1,:],self.elements)   
    
        pyplot.figure()            # plot triangulation
        pyplot.gca().set_aspect('equal')
        pyplot.triplot(triang,'k-',lw=1)
        pyplot.title('Some mesh', fontsize = 15)  
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)
        pyplot.show()

    def plotz(self,zfun):          # plots the level-set function as surface and as contours
        triang = tri.Triangulation(self.points[0,:], self.points[1,:],self.elements)
        z      = zfun(self.points[0,:],self.points[1,:])
        
        pyplot.figure()            # plots countours of some function zfun
        pyplot.gca().set_aspect('equal')
        pyplot.tricontourf(triang, z)
        pyplot.colorbar()
        pyplot.tricontour(triang,z,colors='k',levels=[-0.5,0,0.5])
        pyplot.title('Some contours', fontsize = 15)  
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)

        fig = pyplot.figure()      # plot surface
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(triang,z, cmap=pyplot.cm.CMRmap,alpha=0.75)
        pyplot.title('Some function', fontsize = 15)  
        ax.set_xlabel('$x$', fontsize = 20)
        ax.set_ylabel('$y$', fontsize = 20)
        ax.set_zlabel('$z$',fontsize = 20)
        
        pyplot.show()              # shows and saves plots
        fig.savefig('Trisurf_exact.svg', dpi=fig.dpi, bbox_inches = "tight")


    def findzero(self,z):          # finds the location of the level-set curve:
        zeroelements   = []        # the triangles of the mesh...
        zeroedges      = []        # the edges of the mesh...
        zeroedgepoints = []        # the points of the edges of the mesh...which it intersects
        zeropoints     = []        # the zero points of the level-set along the zeroedges
        levelsetedges  = []        # the level-set curve through the zeroelement
        
        for e in self.edges:                         # finds zeroedges
            x1 = self.points[0,e[0]]
            y1 = self.points[1,e[0]]
            x2 = self.points[0,e[1]]
            y2 = self.points[1,e[1]]
            if z(x1,y1)*z(x2,y2)<0:
                zeroedges.append(e)
##                for p in [e[0],e[1]]:              # lists the extremes of zeroedges
##                    if zeroedgepoints == []:
##                        zeroedgepoints.append(p)
##                    else:
##                        i=0
##                        while p < zeroedgepoints[i]:
##                            i+=1
##                        if not p == zeroedgepoints[i]:
##                            zeroedgepoints.insert(i+1,p)
                if x1 == x2:                         # finds the zeropoints along zeroedges
                    x0 = x1
                else:
                    x0 = x1+(x2-x1)*z(x1,y1)/(z(x1,y1)-z(x2,y2))
                    if x0 > x1 and x0 > x2:
                        x0 = max(x1,x2)
                    elif x0 < x1 and x0 < x2:
                        x0 = min(x1,x2)
                if y1 == y2:
                    y0 = y1
                else:
                    y0 = y1+(y2-y1)*z(x1,y1)/(z(x1,y1)-z(x2,y2))
                    if y0 > y1 and y0 > y2:
                        y0 = max(y1,y2)
                    elif y0 < y1 and y0 < y2:
                        y0 = min(y1,y2)
                zeropoints.append([x0,y0])

        for el in self.elements:                    # finds the zeroelements
            n = 0
            for e in zeroedges:
                if e[0] in el and e[1] in el:
                    n+=1
            if n == 2:
                zeroelements.append(el)
            elif n == 1:
                z1 = z(self.points[0,el[0]],self.points[1,el[0]])
                z2 = z(self.points[0,el[1]],self.points[1,el[1]])
                z3 = z(self.points[0,el[2]],self.points[1,el[2]])
                if z1*z2*z3 == 0:
                    zeroelements.append(el)
            if n > 2:
                print('------------error------------- too many zeroedges around element')
            
        for el in zeroelements:                     # finds the levelsetedges
            p = []
            for i in range(len(zeroedges)):
                if zeroedges[i][0] in el and zeroedges[i][1] in el:
                    p.append(zeropoints[i])
            levelsetedges.append(p)

        self.zeroelements   = zeroelements
        self.zeroedges      = zeroedges
        self.zeroedgepoints = zeroedgepoints
        self.zeropoints     = zeropoints
        self.levelsetedges  = levelsetedges
        return

    def altfindzero(self,z):       # finds the location of the level-set curve:
        self.zeroelements   = []        # the triangles of the mesh...
        self.zeroedges      = []        # the edges of the mesh...which the level-set curve passes through
        tempzeroedges       = []
        self.zeropoints     = []        # the zero points of the level-set along the zeroedges
        self.levelsetedges  = []        # the level-set curve through the zeroelement

        for elem in self.elements:
            p1 = [self.points[0][elem[0]],self.points[1][elem[0]]]
            p2 = [self.points[0][elem[1]],self.points[1][elem[1]]]
            p3 = [self.points[0][elem[2]],self.points[1][elem[2]]]
            pxyarr = [p1,p2,p3]
            numzeroedges = 0
            for pair in [[0,1],[0,2],[1,2]]:
                x1 = pxyarr[pair[0]][0]
                y1 = pxyarr[pair[0]][1]
                x2 = pxyarr[pair[1]][0]
                y2 = pxyarr[pair[1]][1]
                if z(x1,y1)*z(x2,y2)<0:                  # finds the zeroedges
                    tempzeroedges.append([min(elem[pair[0]],elem[pair[1]]),max(elem[pair[0]],elem[pair[1]])])
                    numzeroedges += 1
                    if x1 == x2:                         # finds the zeropoints along zeroedges
                        x0 = x1
                    else:
                        x0 = x1+(x2-x1)*z(x1,y1)/(z(x1,y1)-z(x2,y2))
                        if x0 > x1 and x0 > x2:
                            x0 = max(x1,x2)
                        elif x0 < x1 and x0 < x2:
                            x0 = min(x1,x2)
                    if y1 == y2:
                        y0 = y1
                    else:
                        y0 = y1+(y2-y1)*z(x1,y1)/(z(x1,y1)-z(x2,y2))
                        if y0 > y1 and y0 > y2:
                            y0 = max(y1,y2)
                        elif y0 < y1 and y0 < y2:
                            y0 = min(y1,y2)
                    self.zeropoints.append([x0,y0])
                    
            if numzeroedges == 2:                        # finds the zeroelements
                self.zeroelements.append(elem)
                self.levelsetedges.append([self.zeropoints[-1],self.zeropoints[-2]])

        zeroedgeindices = findedgesindices(tempzeroedges,self.edges)    # finds the indices of the zeroedges
        self.zeroedges = zeroedgeindices

        return

    def plotzeroz(self):                 # plots the mesh along with the zeroelements, edges and points
        triang = tri.Triangulation(self.points[0,:], self.points[1,:],self.elements)  # ...found in findzero
    
        pyplot.figure()

        pyplot.gca().set_aspect('equal') # plot triangulation
        pyplot.triplot(triang,'k-',lw=1)
        pyplot.title('Some mesh', fontsize = 15)  
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)

        for el in self.zeroelements:     # fill the zeroelements light grey
            x = [self.points[0,el[0]],self.points[0,el[1]],self.points[0,el[2]]]
            y = [self.points[1,el[0]],self.points[1,el[1]],self.points[1,el[2]]]
            pyplot.fill(x,y,'xkcd:light grey')
        for i in self.zeroedges:         # plot the zeroedges red
            p1 = self.edges[i][0]
            p2 = self.edges[i][1]
            x = [self.points[0,p1],self.points[0,p2]]
            y = [self.points[1,p1],self.points[1,p2]]
            pyplot.plot(x,y,'r-',lw=1.2)
        for e in self.levelsetedges:     # plot the levelsetedges blue
            x = [e[0][0],e[1][0]]
            y = [e[0][1],e[1][1]]
            pyplot.plot(x,y,'xkcd:bright blue')
        for p in self.zeropoints:        # plot the zeropoints black
            pyplot.plot(p[0],p[1],'ko',markersize = 1.5)
        pyplot.show()

    
        
                

