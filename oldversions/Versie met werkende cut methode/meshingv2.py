# -*- coding: utf-8 -*-
"""
Created on 

@author: twortelboer, douden
"""

import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.tri as tri

import matlab
import matlab.engine as me
from mpl_toolkits.mplot3d import Axes3D



#### SOME FUNCTIONS ############################################################################

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

def pythagoras(p1,p2):
    if isinstance(p1,point):
        p1coords = p1.coordinates
    else:
        p1coords = p1
    if isinstance(p2,point):
        p2coords = p2.coordinates
    else:
        p2coords = p2
    return np.sqrt( (p1coords[0]-p2coords[0])**2 + (p1coords[1]-p2coords[1])**2 )

def zzerointerp(p1,p2,zfun):  # finds the zeropoints along zeroedges
    [x1,y1] = p1
    [x2,y2] = p2
    zp1 = zfun([x1,y1])
    zp2 = zfun([x2,y2])
    if x1 == x2:
        x0 = x1
    else:
        x0 = x1+(x2-x1)*abs(zp1/(zp1-zp2))
    if y1 == y2:
        y0 = y1
    else:
        y0 = y1+(y2-y1)*abs(zp1/(zp1-zp2))
    return [x0,y0]

def edgescommonp(e1,e2):
    for p in e1.points:
        if p in e2.points:
            return p

def edgeswithpoint(p,elist):
    edges = []
    for e in elist:
        if p in e:
            edges.append(e)
    return edges


#### INITIALIZATION OF point class ################################################################

class point:                           # Builds a point class consisting of coordinates
    def __init__(self,x,y,mesh):
        self.coordinates = [x,y]
        self.nbedges     = []
        self.nbelements  = []
        self.iszeropoint = False
        mesh.points.append(self)

    def __repr__(self):                 # Adds a string representation with some basic info
        string  = 'Point object with \n'
        string += 'Coordinates: x = ' + str(self.coordinates[0]) + '\n'
        string += '             y = ' + str(self.coordinates[1]) + '\n'
        string += 'Number of incident edges = ' + str(len(self.nbedges)) + '\n'
        string += 'Number of indicent elements = ' + str(len(self.nbedges))
        return string

    def __eq__(self,other):
        if self is other:
            return True
        else:
            return False

    def delete(self,mesh):
        for e in self.nbedges:
            e.delete(mesh)
        mesh.points = [p for p in mesh.points if p is not self]
        del self
        return

    def move(self,p,z):     # Performs the necessary operations to move a point to new coordinates p
        self.coordinates[0] = p[0]
        self.coordinates[1] = p[1]

        for e in self.nbedges:
            e.length = e.calclength()
            e.findzero(z)
        for el in self.nbelements:
            el.size = el.calcsize()
            el.skewness = el.calcskew()
            el.findzero(z)
        return


#### INITIALIZATION OF edge class ################################################################

class edge:                        # Builds an edge class consisting of two points
    def __init__(self,p1,p2,mesh):
        self.points     = ( p1, p2 )
        for p in self.points:
            p.nbedges += (self,)
        self.nbelements = []
        self.length     = self.calclength()
        self.islevelset = False
        self.iszeroedge = False
        mesh.edges.append(self)

    def __repr__(self):
        string  = 'Edge object with \n'
        string += 'Number of incident elements = ' + str(len(self.nbelements)) + '\n'
        return string

    def __eq__(self,other):              # Decides if some other line is equal based on point objects
        (s1,s2) = self.points
        if isinstance(other,edge):
            if self is other:
                return True
            else:
                return False
        else:
            (o1,o2) = other
            if ( (s1 == o1) and (s2 == o2) ) or ( (s1 == o2) and (s2 == o1) ):
                return True
            else:
                return False

    def delete(self,mesh,delobject = True):
        for el in self.nbelements:
            el.delete(mesh)
        for p in self.points:
            p.nbedges = [e for e in p.nbedges if e is not self]
        mesh.edges = [e for e in mesh.edges if e is not self]
        if delobject:
            del self
        return
            

    def calclength(self):                # Calculates the length of the line
        return pythagoras( self.points[0], self.points[1] )

    def findzero(self,zfun):             # Finds out whether the zero-levelset line crosses the edge
                                         # ... and returns the interpolated position of the crossing
        p1coords = self.points[0].coordinates
        p2coords = self.points[1].coordinates
        fp1 = zfun( p1coords )
        fp2 = zfun( p2coords )
        testval = fp1 * fp2
        if testval < 0 and abs(fp1) > self.hstep/10 and abs(fp2) > self.hstep/10:
            self.iszeroedge = True
            self.zeropoint  = zzerointerp( p1coords, p2coords, zfun )
            self.islevelset = False
            return
        else:
            self.iszeroedge = False
            n=0
            if abs(fp1) < self.hstep/10:
                self.points[0].iszeropoint = True
                n += 1
            else:
                self.points[0].iszeropoint = False
            if abs(fp2) < self.hstep/10:
                self.points[1].iszeropoint = True
                n += 1
            else:
                self.points[1].iszeropoint = False
            if n == 2:
                self.islevelset = True
            else:
                self.islevelset = False
            return
        
    def otherp(self,selfp):
        for p in self.points:
            if p is selfp:
                None
            else:
                return p

    def split(self,mesh,deledgeobject = True):              # Creates a new point on the zeropoint and two new lines
        p1 = self.points[0]            # ... and returns the created elements.
        p2 = self.points[1]
        zp = self.zeropoint

        np  = point( zp[0], zp[1], mesh )
        ne1 = edge( np, p1, mesh )
        ne2 = edge( np, p2, mesh )

        self.delete(mesh,delobject = deledgeobject)
        
        return [np,ne1,ne2]




#### INITIALIZATION OF element class ################################################################

class element:                           # Builds an element class consisting of three points and three lines
    def __init__(self, p1, p2, p3, mesh, edges = None):
        self.points = ( p1, p2, p3 )
        for p in self.points:
            p.nbelements.append(self)

        if not edges == None:
            self.edges = edges
        else:
            self.edges  = []
            for e in p1.nbedges:
                if e == [p1,p2] or e == [p1,p3]:
                    self.edges.append(e)
                    e.nbelements.append(self)
            for e in p2.nbedges:
                if e == [p2,p3]:
                    self.edges.append(e)
                    e.nbelements.append(self)
        self.edges = tuple(self.edges)
        self.nbelements = []
        for e in self.edges:
            for el in e.nbelements:
                if el is not self:
                    self.nbelements.append(el)
        
        if len(self.edges)<3:
            print('Error not enough edges found in init method of element')
            for p in self.points:
                print(p)
            for e in self.edges:
                print(e)

        self.size     = self.calcsize()
        self.skewness = self.calcskew()
        mesh.elements.append(self)


    def __repr__(self):
        (p1,p2,p3) = self.points
        string  = 'Element object with \n'
        string += 'The coordinates of point 1 are: ' + str(p1.coordinates) + '\n'
        string += '        ,,         point 2 are: ' + str(p2.coordinates) + '\n'
        string += '        ,,         point 3 are: ' + str(p3.coordinates) + '\n'
        n = 1
        if len(self.zpointonedge) > 0:
            substring =     'The coordinates of zpoint ' + str(n) + ' on an edge are: '
            for zp in self.zpointonedge:
                string += substring + str(zp) + '\n'
                n += 1
                substring = '         ,,        zpoint ' + str(n) + ' on an edge are: '
        if len(self.zpoints) > 0:
            substring =     'The coordinates of zpoint ' + str(n) + ' on an edge are: '
            for zp in self.zpoints:
                string += substring + str(zp.coordinates)
                n += 1
                substring = '         ,,        zpoint ' + str(n) + ' on an edge are: '
        string += 'Points:\n'
        for p in self.points:
            string += repr(p)
        string += '\n Edges:\n'
        for e in self.edges:
            string += repr(e)
        return string

    def __eq__(self,other):
        if isinstance(other,element):
            if self is other:
                return True
            else:
                return False

    def delete(self,mesh, delobject = True):
        for el in self.nbelements:
            try:
                el.nbelements.remove(self)
            except:
                None
        for p in self.points:
            p.nbelements = [el for el in p.nbelements if el is not self]
        for e in self.edges:
            e.nbelements = [el for el in e.nbelements if el is not self]
        mesh.elements = [el for el in mesh.elements if el is not self]
        if delobject:
            del self
        return

    def findmid(self):
        (p1,p2,p3) = self.points
        (x1,y1) = p1.coordinates
        (x2,y2) = p2.coordinates
        (x3,y3) = p3.coordinates
        return ( (x1+x2+x3)/3, (y1+y2+y3)/3 )

    def calcskew(self):                        # Calculates the skewness of the element
        d1 = self.edges[0].length
        d2 = self.edges[1].length
        d3 = self.edges[2].length
        a1 = np.arccos( (d1**2+d2**2-d3**2) / (2*d1*d2) )
        a2 = np.arccos( (d2**2+d3**2-d1**2) / (2*d2*d3) )
        a3 = np.pi - a1 - a2
        amax = max(a1,a2,a3)
        amin = min(a1,a2,a3)
        return max( (amax - np.pi/3) / (2*np.pi/3) , (np.pi/3 - amin) / (np.pi/3) )

    def calcsize(self):                         # Calculates the size of the element
        d1 = self.edges[0].length
        d2 = self.edges[1].length
        d3 = self.edges[2].length
        p  = (d1+d2+d3)/2
        return np.sqrt( p*(p-d1)*(p-d2)*(p-d3) )

    def otherp(self,p1,p2):
        for p in self.points:
            if (p == p1) or (p == p2):
                None
            else:
                return p

    def findzero(self,zfun):
        self.zedges       = ()                  # Edges crossing the zero-levelset line
        self.zpointonedge = ()                  # Interpolated zero-levelset point on edge
        self.levelsetlines = []
        
        nze               = 0                   # Look for crossings of zero-levelset line with edges
        for e in self.edges:
            if e.iszeroedge == True:
                self.zedges       += (e,)
                nze               += 1
                self.zpointonedge += (e.zeropoint,)

        nzp          = 0                        # Look for element points on the zero-levelset line
        self.zpoints = ()
        for p in self.points:
            if p.iszeropoint == True:
                self.zpoints      += (p,)
                nzp               += 1

        if nze == 0 and nzp == 0:               # Return false when no zero-levelset line incident
            self.iszeroelement = False
            return
        elif nzp == 0:
            if nze == 2:                        # ... true when line crosses two edges
                self.iszeroelement = True
                self.levelsetlines = [ list(self.zpointonedge) ]
                return
            elif nze == 3:                      # ... true when line crosses three edges
                self.iszeroelement = True
                self.midpoint = self.findmid()
                self.levelsetlines = [ (self.zpointonedge[0],self.midpoint),
                                       (self.zpointonedge[1],self.midpoint),
                                       (self.zpointonedge[2],self.midpoint) ]
                return
            else:
                raise NameError('in findzero of element; nzp=0 and not nze=2 or 3')
        elif nzp == 1:
            if nze == 1:                        # ... true when line crosses an edge and a line
                self.iszeroelement = True
                self.levelsetlines = [ (self.zpointonedge[0], self.zpoints[0].coordinates) ]
                return 
            elif nze == 0:                      # ... false when line incident on just one point
                self.iszeroelement = False
                return
            elif nze == 2:
                self.iszeroelement = True
                self.midpoint = self.findmid()
                self.levelsetlines = [ (self.zpointonedge[0],self.midpoint), (self.zpointonedge[1],self.midpoint), (self.zpoints[0].coordinates,self.midpoint) ]
                return
            else:
                raise NameError('in findzero of element; nzp=1 and not nze=1 or 2')
        elif nzp == 2:
            if nze == 0:                        # ... false when line incident on just two points
                self.iszeroelement = False
                return
            elif nze == 1:                      # ... true when line crosses an edge and two points
                self.iszeroelement = True
                self.midpoint = self.findmid()
                self.levelsetlines = [ (self.zpointonedge[0],self.midpoint),
                                       (self.zpoints[0].coordinates,self.midpoint),
                                       (self.zpoints[1].coordinates,self.midpoint) ]
                return
            elif nze == 2:
                self.iszeroelement = True
                self.midpoint = self.findmid()
                self.levelsetlines = [ (self.zpointonedge[0],self.midpoint),
                                       (self.zpointonedge[1],self.midpoint),
                                      (self.zpoints[0].coordinates,self.midpoint),
                                       (self.zpoints[1].coordinates,self.midpoint) ]
                return
            else:
                raise NameError('in findzero of element; nzp=2 and not nze=0 or 1')
        elif nzp == 3:                          # ... true when line incident on all three points
            self.iszeroelement = True
            self.midpoint = self.findmid()
            (p1,p2,p3) = self.points
            self.levelsetlines = [ (p1.coordinates,self.midpoint),
                                   (p2.coordinates,self.midpoint),
                                   (p3.coordinates,self.midpoint) ]
            return
        else:
            raise NameError('in findzero of element, not nzp=0,1 or 2')

    def cut( self, mesh ):
        newcomponents = []
        splitedges = []
        
        self.delete(mesh,delobject = False)
        
        for e in self.edges:
            if e.iszeroedge == True:
                newcomponents.append( e.split(mesh,deledgeobject = False) )
                splitedges.append(e)
                
        if len(newcomponents) == 1:
            [np,ne1,ne2] = newcomponents[0]
            (p1,p2) = splitedges[0].points
            op = self.otherp(p1,p2)
            otheredges = list(self.edges[:])
            for e in otheredges:
                if e == (p1,p2):
                    otheredges.remove(e)
                    break
            if otheredges[0] == (op,p1):
                oe1 = otheredges[0]
                oe2 = otheredges[1]
            else:
                oe1 = otheredges[1]
                oe2 = otheredges[0]

            ne = edge( op, np, mesh )
            ne.islevelset = True
            element( op, np, p1, mesh, edges = (ne,ne1,oe1) )
            element( op, np, p2, mesh, edges = (ne,ne2,oe2) )
            
        if len(newcomponents) == 2:
            if splitedges[0].points[0] in splitedges[1].points:
                cp = splitedges[0].points[0]
            else:
                cp = splitedges[0].points[1]
            if pythagoras( cp, splitedges[0].zeropoint ) > pythagoras( cp, splitedges[1].zeropoint ):
                [sp,se1,se2] = newcomponents[0]
                [np,ne1,ne2] = newcomponents[1]
                op1 = self.otherp( splitedges[0].points[0], splitedges[0].points[1] )
                op2 = self.otherp( splitedges[1].points[0], splitedges[1].points[1] )
            else:
                [sp,se1,se2] = newcomponents[1]
                [np,ne1,ne2] = newcomponents[0]
                op1 = self.otherp( splitedges[1].points[0], splitedges[1].points[1] )
                op2 = self.otherp( splitedges[0].points[0], splitedges[0].points[1] )
            if cp in se2.points:
                dummy = se1
                se1 = se2
                se2 = dummy
            if cp in ne2.points:
                dummy = ne1
                ne1 = ne2
                ne2 = dummy
            for e in self.edges:
                if not cp in e.points:
                    re = e
            de1 = edge( op1, sp, mesh )
            de2 = edge( sp, np, mesh )
            de2.islevelset = True
            el1 = element( sp, op1, op2, mesh, edges = (re,se2,de1) )
            el2 = element( sp, op1, np, mesh, edges = (de1,de2,ne2) )
            el3 = element( sp, np, cp, mesh, edges = (se1,ne1,de2) )

        if len(newcomponents) > 2:
            print('Error --- too many zeroedges for cut method.')

        for e in splitedges:
            del e
        del self
        
        return

#### INITIALIZATION OF mesh class ######################################################################

class mesh:

    def __init__(self,hstep,threshold,
                     x = np.array([-1,1,1,-1]),
                     y = np.array([-1,-1,1,1])  ):
        h = np.array(hstep)
        self.hstep = hstep
        self.threshold = threshold
            # check whether saved mesh exists
        file_str = "h"+np.array2string(h)+"_x"+np.array2string(x)+"_y"+np.array2string(y)+".npz"
        try:
            npzfile = np.load(file_str)
            for key in npzfile.keys():
                exec("self.f" + key + " = npzfile[key]")
        except IOError:
            print("Mesh does not exist\nCreating the mesh")
            # create mesh, random numbers
            self.felements, self.fedges, self.fpoints = create_mesh(x,y,h)
            print("Saving the mesh")
            np.savez(file_str, elements=elements,edges=edges,points=points)
            print("Mesh saved")

        self.points = []
        for i,x in enumerate(self.fpoints[0]):
            point( x, self.fpoints[1][i], self )
            
        self.edges = []
        for e in self.fedges:
            newedge = edge( self.points[e[0]], self.points[e[1]], self )
            newedge.threshold = threshold
            newedge.hstep = hstep

        self.elements = []
        for el in self.felements:
            p1 = self.points[el[0]]
            p2 = self.points[el[1]]
            p3 = self.points[el[2]]
            newelem = element( p1, p2, p3, self )
            newelem.threshold = threshold
            newelem.hstep = hstep

    def initzeroedgesinfo(self,zfun):
        for e in self.edges:
            e.findzero(zfun)
        return

    def findzero(self,zfun):       # lists zeroelements for plotting purposes
        self.levelsetlines = []    # lists the zero-levelset lines not on the grid
        for el in self.elements:
            if el.iszeroelement == True:
                self.levelsetlines.extend( el.levelsetlines )

        self.zeropoints = ()
        for p in self.points:
            if p.iszeropoint == True:
                self.zeropoints += (p,)
        return

    def updatezeroelementinfo(self,zfun):
        for el in self.elements:
            el.findzero(zfun)
        return


#### MESH ADJUSTMENT METHODS ##########################################################################

    def simplegridtozeropoint(self,z,thresh):
        self.zeroedges = ()
        for e in self.edges:
            if e.iszeroedge == True:
                self.zeroedges += (e,)

        adjustments = []
        for ze in self.zeroedges:
            p1      = ze.points[0]
            p2      = ze.points[1]
            (x1,y1) = p1.coordinates
            (x2,y2) = p2.coordinates
            zp      = ze.zeropoint
            p1dist  = pythagoras(p1,zp)
            if p1dist < thresh:
                if x1 > -1 and x1 < 1 and y1 > -1 and y1 < 1:
                    adjustments.append([p1,zp,p1dist])
                elif (x1 == -1 or x1 == 1) and zp[0] == x1:
                    adjustments.append([p1,zp,p1dist])
                elif (y1 == -1 or y1 == 1) and zp[1] == y1:
                    adjustments.append([p1,zp,p1dist])
            p2dist = pythagoras(p2,zp)
            if p2dist < thresh:
                if x2 > -1 and x2 < 1 and y2 > -1 and y2 < 1:
                    adjustments.append([p2,zp,p2dist])
                elif x2 == -1 or x2 == 1 and zp[0] == x2:
                    adjustments.append([p2,zp,p2dist])
                elif y2 == -1 or y2 == 1 and zp[1] == y2:
                    adjustments.append([p2,zp,p2dist])
        adjustments.sort(key=lambda x: id(x[0]), reverse=True)      # making use of the sorted adjustments
        if len(adjustments)>0:                                      # list to choose most favorable
            adjustments2 = [adjustments[0]]                         # adjustment based on distance
            for i in range(1,len(adjustments)):
                if not adjustments[i][0] is adjustments2[-1][0]:
                    adjustments2.append(adjustments[i])
                elif adjustments[i][2] < adjustments[i-1][2]:
                    adjustments2[-1] = adjustments[i]
        else:
            adjustments2 = []

        for adj in adjustments2:                                    # making the adjustment
            adj[0].move(adj[1],z)
        return

    def newtriang(self):
        for el in self.elements[:]:
            if el.iszeroelement == True:
                el.cut(self)
        return


#### QUALITY METHODS ###################################################################################

    def skewness(self):
        skewness = []
        for el in self.elements:
            skewness.append( el.skewness )
        return sum(skewness)/len(skewness), max(skewness), np.std(skewness)

    def sizes(self):
        sizes = []
        for el in self.elements:
            sizes.append( el.size )
        return min(sizes), max(sizes), np.std(sizes)


#### PLOT METHODS #######################################################################################

    def converttocoordlist(self):
        xcoords = []
        ycoords = []
        for i,p in enumerate(self.points):
            xcoords.append( p.coordinates[0] )
            ycoords.append( p.coordinates[1] )
            p.index = i

        elemindiceslist = []
        for el in self.elements:
            indices = []
            for p in el.points:
                indices.append( p.index )
            elemindiceslist.append( indices )

        return ( [xcoords, ycoords], elemindiceslist )

    def plot(self, show = True):
        (points, elements) = self.converttocoordlist()
        triang = tri.Triangulation(points[0], points[1], elements)

        pyplot.figure()            # plot triangulation
        pyplot.gca().set_aspect('equal')
        pyplot.triplot(triang,'k-',lw=1) 
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)
        if show:
            pyplot.show()

    def plotz(self, zfun):    # plots the level-set function as surface and as contours
        (points, elements) = self.converttocoordlist()
        triang = tri.Triangulation(points[0], points[1], elements)

        z = zfun([points[0],points[1]])
        z = np.array(z)
        
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


    def plotzeroz(self, savestring):                 # plots the mesh along with the zeroelements, edges and points
##        (points, elements) = self.converttocoordlist()
##        triang = tri.Triangulation(points[0], points[1], elements) # ...found in findzero

        fig = pyplot.figure(figsize=(10,10))
        pyplot.gca().set_aspect('equal') # plot triangulation
        #pyplot.triplot(triang,'k-',lw=1)  
        pyplot.xlabel('$x$', fontsize = 20)
        pyplot.ylabel('$y$', fontsize = 20)
        pyplot.xticks(fontsize = 15)
        pyplot.yticks(fontsize = 15)

        for el in self.elements:     # fill the zeroelements light grey
            if el.iszeroelement == True:
                (p1,p2,p3) = el.points
                (x1,y1) = p1.coordinates
                (x2,y2) = p2.coordinates
                (x3,y3) = p3.coordinates
                x = [x1,x2,x3]
                y = [y1,y2,y3]
                pyplot.fill(x,y,'xkcd:light grey')

        for e in self.edges:            # plots the edges according to type
            (p1,p2) = e.points
            (x1,y1) = p1.coordinates
            (x2,y2) = p2.coordinates
            x       = [x1,x2]
            y       = [y1,y2]
            if e.islevelset == True:              
                pyplot.plot(x,y,'xkcd:bright green')
            elif e.iszeroedge == True:
                pyplot.plot(x,y,'r-',lw=1.2)
            else:
                pyplot.plot(x,y,'k-',lw=1)

        for e in self.levelsetlines:     # plot the levelsetedges blue
            (p1,p2) = ( e[0], e[1] )
            (x1,y1) = p1
            (x2,y2) = p2
            x       = [x1,x2]
            y       = [y1,y2]
            pyplot.plot(x,y,'xkcd:bright blue')

        for p in self.zeropoints:        # plot the zeropoints black
            (x,y) = p.coordinates
            pyplot.plot(x,y,'ko',markersize = 1.5)

        fig.savefig(savestring, dpi=fig.dpi, bbox_inches = "tight")

            



