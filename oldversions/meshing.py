# -*- coding: utf-8 -*-
"""
Created on 

@author: twortelboer, douden
"""

#from mpl_toolkits.mplot3d import Axes3D
import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.tri as tri
##matplotlib.rcParams['text.usetex'] = True
##from matplotlib import rc
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#### for Palatino and other serif fonts use:
###rc('font',**{'family':'serif','serif':['Palatino']})
##rc('text', usetex=True)

import matlab
import matlab.engine as me
from mpl_toolkits.mplot3d import Axes3D
#from copy import deepcopy as dcopy

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

def findedgesindices(edges, edgelist):  # finds the indices of edges in edgelist
    edges.sort()                        # based on point indices
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
        if x0 > x1 and x0 > x2:          # makes sure x0 is still in (x1,x2)
            x0 = max(x1,x2)
        elif x0 < x1 and x0 < x2:
            x0 = min(x1,x2)
    if y1 == y2:
        y0 = y1
    else:
        y0 = y1+(y2-y1)*z(x1,y1)/(z(x1,y1)-z(x2,y2))
        if y0 > y1 and y0 > y2:          # makes sure y0 is still in (y1,y2)
            y0 = max(y1,y2)
        elif y0 < y1 and y0 < y2:
            y0 = min(y1,y2)
    return [x0,y0]

def edgeisequal(e1,e2):
    if e1[0] == e2[0] and e1[1] == e2[1]:
        return True
    else:
        return False

def elemisequal(elem1,elem2):
    if elem1[0] in elem2 and elem1[1] in elem2 and elem1[2] in elem2:
        return True
    else:
        return False

def edgeinelem(e,elem,findopp = False):
    if findopp:
        for i in range(3):
            if not elem[i] in e:
                return [True,elem[i]]
    else:    
        if e[0] in elem and e[1] in elem:
            return True
        else:
            return False

class meshclass:
    def __init__(self, hstep,
                     x = np.array([-1,1,1,-1]),
                     y = np.array([-1,-1,1,1])   ):
        h = np.array(hstep)
        self.hstep = hstep
            # check whether saved mesh exists
        file_str = "h"+np.array2string(h)+"_x"+np.array2string(x)+"_y"+np.array2string(y)+".npz"
        try:
            npzfile = np.load(file_str)
            for key in npzfile.keys():
                exec("self."+ key + " = npzfile[key]")
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
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)
        pyplot.show()

    def plotz(self,zfun):          # plots the level-set function as surface and as contours
        triang = tri.Triangulation(self.points[0,:], self.points[1,:],self.elements)
        z      = zfun(self.points[0,:],self.points[1,:])
        z = np.array(z)
        print(str(len(z)))
        
        fig = pyplot.figure(figsize=(10,10))            # plots countours of some function zfun
        pyplot.gca().set_aspect('equal')
        pyplot.tricontourf(triang, z)
        pyplot.colorbar()
        pyplot.tricontour(triang,z,colors='k',levels=[-0.5,0,0.5]) 
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)
        
        string = 'Saved_figures/' + zfun.__name__ + '_contour.png'
        fig.savefig(string, dpi=fig.dpi, bbox_inches = "tight")

        fig = pyplot.figure(figsize=(10,10))      # plot surface
        ax = fig.gca(projection='3d')
        ax.plot_trisurf(triang,z, cmap=pyplot.cm.CMRmap,alpha=0.75) 
        ax.set_xlabel('$x$', fontsize = 20)
        ax.set_ylabel('$y$', fontsize = 20)
        ax.set_zlabel('$z$',fontsize = 20)
        
        string = 'Saved_figures/' + zfun.__name__ + '_3dsurf.png'
        fig.savefig(string, dpi=fig.dpi, bbox_inches = "tight")

    def findzero(self,z):       # finds the location of the level-set curve:
        self.zeroelements    = []        # the triangles of the mesh...
        self.zeroedges       = []        # the edges of the mesh...which the level-set curve passes through
    #    self.zeroedgeindices = []
        self.zeroedges       = []
        self.zeropoints      = []        # the zero points of the level-set along the zeroedges
        self.levelsetedges   = []        # the level-set curve through the zeroelement
        zpointongrid         = []
        levelsetongrid       = []
        thresh = 0.01*self.hstep

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
                
                if z(x1,y1)*z(x2,y2) < 0 and abs(z(x1,y1))>thresh and abs(z(x2,y2))>thresh:
                    numzeroedges += 1                                          # finds the zeroedges
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
                                    break
                            if not (self.zeroedges[i][1] == q2 or i > len(self.zeroedges)-1):
                                self.zeroedges.insert(i,[q1,q2])
                                self.zeropoints.insert(i,zpoint[-1])

            if numzeroedges == 2:                        # finds the zeroelements
                self.zeroelements.append(elem)
                self.levelsetedges.append([zpoint[0],zpoint[1]])
            elif numzeroedges == 1:
                for p in [p1,p2,p3]:
                    if not abs(z(p[0],p[1])) > thresh:
                        zpointongrid.append(p)
                        levelsetongrid.append([p,zpoint[-1]])
                        self.zeroelements.append(elem)
                        break
            elif numzeroedges == 0:
                zp = []
                for p in [p1,p2,p3]:
                    if not abs(z(p[0],p[1])) > thresh:
                        zp.append(p)
                if len(zp) == 2:
                    zpointongrid.append(zp[0])
                    zpointongrid.append(zp[1])
                    levelsetongrid.append([zp[0],zp[1]])
                    self.zeroelements.append(elem)
                        
                    
        self.zeropoints.extend(zpointongrid)
        self.levelsetedges.extend(levelsetongrid)
 #       self.zeroedgeindices = findedgesindices(self.zeroedges,self.edges)    # finds the indices of the zeroedges   STILL NECECARRY???
        return

    def plotzeroz(self, savestring):                 # plots the mesh along with the zeroelements, edges and points
        triang = tri.Triangulation(self.points[0,:], self.points[1,:],self.elements)  # ...found in findzero

        fig = pyplot.figure(figsize=(10,10))
        pyplot.gca().set_aspect('equal') # plot triangulation
        pyplot.triplot(triang,'k-',lw=1)  
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
        fig.savefig(savestring, dpi=fig.dpi, bbox_inches = "tight")

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
                if x1 > -1 and x1 < 1 and y1 > -1 and y1 < 1:
                    adjustments.append([p1,zp,p1dist])
                elif (x1 == -1 or x1 == 1) and zp[0] == x1:
                    adjustments.append([p1,zp,p1dist])
                elif (y1 == -1 or y1 == 1) and zp[1] == y1:
                    adjustments.append([p1,zp,p1dist])
            p2dist = np.sqrt((x2-zerox)**2 + (y2-zeroy)**2)         
            if p2dist < thresh:
                if x2 > -1 and x2 < 1 and y2 > -1 and y2 < 1:
                    adjustments.append([p2,zp,p2dist])
                elif x2 == -1 or x2 == 1 and zp[0] == x2:
                    adjustments.append([p2,zp,p2dist])
                elif y2 == -1 or y2 == 1 and zp[1] == y2:
                    adjustments.append([p2,zp,p2dist])
        adjustments.sort()                                          # making use of the sorted adjustments
        if len(adjustments)>0:                                      # list to choose most favorable
            adjustments2 = [adjustments[0]]                         # adjustment based on distance
            for i in range(1,len(adjustments)):
                if not adjustments[i][0] == adjustments2[-1][0]:
                    adjustments2.append(adjustments[i])
                elif adjustments[i][2] < adjustments[i-1][2]:
                        adjustments2[-1] = adjustments[i]
        else:
            adjustments2 = []

        for adj in adjustments2:                                    # making the adjustment
            self.points[0][adj[0]] = adj[1][0]
            self.points[1][adj[0]] = adj[1][1]
        return
    
    def newtriang(self,thresh):
        ptoadd = []
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
            p2dist = np.sqrt((x2-zerox)**2 + (y2-zeroy)**2)
            if ( not p1dist < thresh ) and ( not p2dist < thresh ):
                ptoadd.append([zp,i])
            
        for i in range(len(ptoadd)):
            print(str(i)+' of the '+str(len(ptoadd)))
            p = ptoadd[i][0]
            pset1 = np.append( self.points[0], p[0] )
            pset2 = np.append( self.points[1], p[1] )
            self.points = [pset1,pset2]
            
            newpindex = len(self.points[0])-1
            
            zedge = self.zeroedges[ ptoadd[i][1] ]
            for elem in self.zeroelements:
                [edgeinel,oppositep] = edgeinelem(zedge,elem,findopp = True)
                if edgeinel:
                    j = 0
                    for el in self.elements:
                        if elemisequal(el, elem):
                            self.elements = np.delete(self.elements,j,0)
                        j += 1

                    self.elements = np.concatenate(( self.elements,np.array([[zedge[0],newpindex,oppositep]]) ))
                    self.elements = np.concatenate(( self.elements,np.array([[zedge[1],newpindex,oppositep]]) ))
                    p1 = min(oppositep,newpindex)
                    p2 = max(oppositep,newpindex)
                    k = 0
                    for e in self.edges:
                        if edgeisequal(zedge,e):
                            self.edges = np.delete( self.edges, k, 0 )
                        k += 1
                    
                    self.edges = np.concatenate((self.edges,np.array([[p1,p2]])))
        self.points = np.array(self.points)
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
        return sum(skewness)/len(skewness), max(skewness), np.std(skewness)

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
