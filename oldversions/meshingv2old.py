# -*- coding: utf-8 -*-
"""
Created on 

@author: twortelboer, douden
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

def zzerointerp(x1,y1,x2,y2,z):
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
    return [x0,y0]

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

    def findzero(self,z):       # finds the location of the level-set curve:
        self.zeroelements    = []        # the triangles of the mesh...
        self.zeroedges       = []        # the edges of the mesh...which the level-set curve passes through
        self.zeroedgeindices = []
        self.zeroedges       = []
        self.zeropoints      = []        # the zero points of the level-set along the zeroedges
        self.levelsetedges   = []        # the level-set curve through the zeroelement
        zpointongrid = []
        levelsetongrid = []

        for elem in self.elements:
            p1 = [self.points[0][elem[0]],self.points[1][elem[0]]]
            p2 = [self.points[0][elem[1]],self.points[1][elem[1]]]
            p3 = [self.points[0][elem[2]],self.points[1][elem[2]]]
            pxyarr = [p1,p2,p3]
            numzeroedges = 0
            zedge = []
            zpoint = []

            for pair in [[0,1],[0,2],[1,2]]:
                x1 = pxyarr[pair[0]][0]
                y1 = pxyarr[pair[0]][1]
                x2 = pxyarr[pair[1]][0]
                y2 = pxyarr[pair[1]][1]

                if z(x1,y1)*z(x2,y2)<0:                  # finds the zeroedges
                    numzeroedges += 1
                    q1 = min(elem[pair[0]],elem[pair[1]])
                    q2 = max(elem[pair[0]],elem[pair[1]])
                    zedge.append([q1,q2])
                    zpoint.append(zzerointerp(x1,y1,x2,y2,z))
                    if self.zeroedges == []:
                        self.zeroedges.append(zedge[-1])
                        self.zeropoints.append(zpoint[-1])
                    else:
                        i = 0
                        while self.zeroedges[i][0] < q1:
                            i += 1
                            if i > len(self.zeroedges)-1:
                                self.zeroedges.append([q1,q2])
                                self.zeropoints.append(zpoint[-1])
                                break
                        if self.zeroedges[i][0] > q1 and i < len(self.zeroedges):
                            self.zeroedges.insert(i,[q1,q2])
                            self.zeropoints.insert(i,zpoint[-1])
                        else:
                            while self.zeroedges[i][0] == q1 and self.zeroedges[i][1] < q2:
                                i += 1
                                if i > len(self.zeroedges)-1:
                                    self.zeroedges.append([q1,q2])
                                    self.zeropoints.append(zpoint[-1])
                            if not (self.zeroedges[i][1] == q2 or i > len(self.zeroedges)-1):
                                self.zeroedges.insert(i,[q1,q2])
                                self.zeropoints.insert(i,zpoint[-1])

            if numzeroedges == 2:                        # finds the zeroelements
                self.zeroelements.append(elem)
                self.levelsetedges.append([zpoint[0],zpoint[1]])
            elif numzeroedges == 1:
                for p in [p1,p2,p3]:
                    if z(p[0],p[1])==0:
                        zpointongrid.append(p)
                        levelsetongrid.append([p,zpoint[-1]])
                        self.zeroelements.append(elem)
                        break
        self.zeropoints.extend(zpointongrid)
        self.levelsetedges.extend(levelsetongrid)
        self.zeroedgeindices = findedgesindices(self.zeroedges,self.edges)    # finds the indices of the zeroedges   STILL NECECARRY???
        return

    def plotzeroz(self,plotshow = True):                 # plots the mesh along with the zeroelements, edges and points
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
        for e in self.zeroedges:         # plot the zeroedges red
            x = [self.points[0,e[0]],self.points[0,e[1]]]
            y = [self.points[1,e[0]],self.points[1,e[1]]]
            pyplot.plot(x,y,'r-',lw=1.2)
        for e in self.levelsetedges:     # plot the levelsetedges blue
            x = [e[0][0],e[1][0]]
            y = [e[0][1],e[1][1]]
            pyplot.plot(x,y,'xkcd:bright blue')
        for p in self.zeropoints:        # plot the zeropoints black
            pyplot.plot(p[0],p[1],'ko',markersize = 1.5)
        if plotshow:
            pyplot.show()

    def simplegridtozeropoint(self,thresh):
        adjustments = []
        for i in range(len(self.zeroedges)):
            p1 = self.zeroedges[i][0]
            p2 = self.zeroedges[i][1]
            x1 = self.points[0][p1]
            y1 = self.points[1][p1]
            x2 = self.points[0][p2]
            y2 = self.points[1][p2]
            zp = self.zeropoints[i]
            zerox = zp[0]
            zeroy = zp[1]
            p1dist = np.sqrt((x1-zerox)**2 + (y1-zeroy)**2)
            if p1dist < thresh:
                adjustments.append([p1,zp,p1dist])
            p2dist = np.sqrt((x2-zerox)**2 + (y2-zeroy)**2)
            if p2dist < thresh:
                adjustments.append([p2,zp,p2dist])
        adjustments = sorted(adjustments)
        adjustments2 = [adjustments[0]]
        for i in range(1,len(adjustments)):
            if not adjustments[i][0] == adjustments2[-1][0]:
                adjustments2.append(adjustments[i])
            elif adjustments[i][2] < adjustments[i-1][2]:
                    adjustments2[-1] = adjustments[i]

        for adj in adjustments2:
            self.points[0][adj[0]] = adj[1][0]
            self.points[1][adj[0]] = adj[1][1]
        return

    def skewness(self):
        skewness = []
        for el in self.elements:
            x1 = self.points[0][el[0]]
            y1 = self.points[1][el[0]]
            x2 = self.points[0][el[1]]
            y2 = self.points[1][el[1]]
            x3 = self.points[0][el[2]]
            y3 = self.points[1][el[2]]
            d1 = np.sqrt( (x2-x1)**2 + (y2-y1)**2 )
            d2 = np.sqrt( (x3-x2)**2 + (y3-y2)**2 )
            d3 = np.sqrt( (x1-x3)**2 + (y1-y3)**2 )
            a1 = np.arccos( (d1**2+d2**2-d3**2) / (2*d1*d2) )
            a2 = np.arccos( (d2**2+d3**2-d1**2) / (2*d2*d3) )
            a3 = np.pi - a1 - a2
            amax = max(a1,a2,a3)
            amin = min(a1,a2,a3)
            skewness.append( max( (amax - np.pi/3) / (2*np.pi/3) , (np.pi/3 - amin) / (np.pi/3) ) )
        return sum(skewness)/len(skewness), max(skewness)

    def sizes(self):
        sizes = []
        for el in self.elements:
            x1 = self.points[0][el[0]]
            y1 = self.points[1][el[0]]
            x2 = self.points[0][el[1]]
            y2 = self.points[1][el[1]]
            x3 = self.points[0][el[2]]
            y3 = self.points[1][el[2]]
            d1 = np.sqrt( (x2-x1)**2 + (y2-y1)**2 )
            d2 = np.sqrt( (x3-x2)**2 + (y3-y2)**2 )
            d3 = np.sqrt( (x1-x3)**2 + (y1-y3)**2 )
            p  = (d1+d2+d3)/2
            sizes.append( np.sqrt(p*(p-d1)*(p-d2)*(p-d3)) )
        return min(sizes),max(sizes),np.std(sizes)
